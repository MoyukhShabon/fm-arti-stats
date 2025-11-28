#!/usr/bin/env python
import polars as pl
import os
import glob
import pysam

## Functions
def return_path_if_exists(path: str, abs=False) -> str:
	if os.path.exists(path):
		return os.path.abspath(path) if abs else path
	else:
		raise FileNotFoundError(f"File not found: {path}")

def read_variants(vcf_path: str, split_multiallelic = True) -> pl.DataFrame:
	variants = (
		pl.read_csv(
			vcf_path, 
			separator="\t", 
			comment_prefix="##", 
			null_values=["."],
			columns=["#CHROM", "POS", "REF", "ALT"], 
			infer_schema_length=100
		).rename(lambda x: x.strip("#").lower())
	)
 
	if split_multiallelic:
		variants = variants.with_columns(
			pl.col("alt").str.split(",").alias("alt")
		).explode("alt")
	
	return variants

def snv_filter(df: pl.DataFrame) -> pl.DataFrame:
	return df.filter(
		(pl.col("ref").str.len_chars() == 1) & 
		(pl.col("alt").str.len_chars() == 1)
	)

def ct_filter(df: pl.DataFrame) -> pl.DataFrame:
	return df.filter(
		((pl.col("ref") == "C") & (pl.col("alt") == "T")) | 
  		((pl.col("ref") == "G") & (pl.col("alt") == "A"))
	)

# Define a function to be applied to each row-like struct
# Used to retrieve neighbourhood sequence of 20bp in each direction via map_elements which is vectorized
def get_neighborhood_sequence(row_struct: dict, genome: pysam.FastaFile) -> str:
    chrom = row_struct["Chr"]
    pos = row_struct["Pos"] # 1-based VCF position
    ref = row_struct["Ref"]
    alt = row_struct["Alt"]

    # Pysam fetch is 0-based, half-open [start, end)
    # VCF Pos 100 = Index 99
    
    # Left Flank (20bp before the variant start)
    # If Pos=100, we want indices 79 to 99 (length 20)
    left_end_idx = pos - 1
    left_start_idx = max(0, left_end_idx - 20)
    
    # Right Flank (20bp after the REF allele ends)
    # If SNV (Ref len 1) at Pos 100: Ref is index 99. Next base is index 100.
    # If Del (Ref len 5) at Pos 100: Ref is 99..103. Next base is index 104.
    right_start_idx = left_end_idx + len(ref)
    right_end_idx = right_start_idx + 20

    try:
        left_flank = genome.fetch(chrom, left_start_idx, left_end_idx)
        right_flank = genome.fetch(chrom, right_start_idx, right_end_idx)

        return (left_flank + alt + right_flank).upper()

    except (ValueError, KeyError):
        # KeyError handles missing chromosomes
        return None

def get_read_length(bam_path: str) -> int:
	bam = pysam.AlignmentFile(bam_path, "rb")
	max_reads = 10000
	
	lengths = []
	for i, read in enumerate(bam):

		if read.is_unmapped:
			continue
		
		if i >= max_reads:
			break
		
		lengths.append(read.query_length)
		
	return pl.Series(lengths).mode()[0]

def natural_sort_variants(df: pl.DataFrame) -> pl.DataFrame:
	return (
		df.with_columns(
			# Create a temporary column 'chr_rank' for sorting
			pl.col("Chr")
			.str.replace("chr", "") # Remove 'chr' prefix
			.str.replace("X", "23") # Handle Sex chromosomes
			.str.replace("Y", "24")
			.str.replace("M", "25") # Handle Mitochondria if present
			.str.replace("MT", "26") # Handle Mitochondria if present
			.cast(pl.Int32, strict=False) # Convert to Integer (strict=False turns unknown contigs to null)
			.fill_null(999) # Put weird contigs at the end
			.alias("chr_rank")
		)
		.sort(["chr_rank", "Pos"]) # Sort by Rank, then Position
		.drop("chr_rank") # Remove the temp column
	)


def annotate_simple_repeats(mut_info: pl.DataFrame, simple_repeats_df: pl.DataFrame) -> pl.DataFrame:
	"""
	Annotates variants with 'Y' or 'N' if they fall within a simple repeat region.
	Assumes both DataFrames contain Chromosome and Position columns.
	"""
	
	# 1. Sort both DataFrames
	# join_asof requires the 'on' column (Pos/start) to be sorted within the 'by' groups (Chr).
	# Sorting by ["Chr", "Pos"] achieves this and keeps the data human-readable.
	variants_sorted = natural_sort_variants(mut_info) #.with_columns(pl.col("Pos").set_sorted())
	repeats_sorted = simple_repeats_df.sort(["chrom", "start"]) #.with_columns(pl.col("start").set_sorted())

	# 2. Perform the ASOF Join
	# We find the repeat segment that started most recently before (or exactly at) the variant position.
	joined = variants_sorted.join_asof(
		repeats_sorted,
		left_on="Pos",
		right_on="start",
		by_left="Chr",
		by_right="chrom",
		strategy="backward" # Look for the closest start position <= variant position
	)

	# 3. Check overlap
	# We found the closest start. Now, is the variant actually INSIDE that segment?
	# i.e., is Variant Pos <= Repeat End?
	final_df = joined.with_columns(
		pl.when(pl.col("Pos") <= pl.col("end"))
		.then(pl.lit("Y"))
		.otherwise(pl.lit("N"))
		.alias("SimpleRepeat_TRF")
	).drop(["start", "end"]) # Clean up BED columns

	return final_df


def make_mut_info(vcf_path: str, genome: pysam.FastaFile, sample_name: str, simple_repeats_bed: pl.DataFrame | None = None, ct_only: bool = False, snv_only: bool = False) -> pl.DataFrame:
	
	# Read Variants
	mut_info = read_variants(vcf_path)

	# Keep only SNVs if requested
	if snv_only:
		mut_info = snv_filter(mut_info)
  
	# Keep only C>T and G>A mutations if requested
	if ct_only:
		mut_info = ct_filter(mut_info)

	# Rename columns to match MicroSEC expected input
	mut_info = mut_info.rename({
		"chrom": "Chr",
		"pos": "Pos",
		"ref": "Ref",
		"alt": "Alt"
	})
 
	# Add Sample Name
	mut_info = mut_info.with_columns(
		pl.lit(sample_name).alias("Sample")
	)

	# MicroSEC requires variants to be at least 200 nucleotides away from the terminal ends of the chromosome/contig where they are located
	chrom_lengths = pl.DataFrame({
		"Chr": genome.references,
		"chr_len_ref": genome.lengths
	})

	# 2. Join and Filter
	mut_info = (
		mut_info.join(chrom_lengths, on="Chr", how="inner")
		.filter(
			(pl.col("Pos") > 200) & 
			(pl.col("Pos") < (pl.col("chr_len_ref") - 200))
		)
		.drop("chr_len_ref")
	)

	# Create a struct of the columns needed for the neighborhood sequence calculation
	# and then apply the function.
	mut_info = mut_info.with_columns(
		pl.struct(["Chr", "Pos", "Ref", "Alt"])
		.map_elements(lambda x: get_neighborhood_sequence(x, genome), return_dtype=pl.Utf8)
		.alias("Neighborhood_sequence")
	)
	
	# Create Mutation Type column as this is a requirement for MicroSEC input
	snv_mask = mut_info["Ref"].str.len_chars() == mut_info["Alt"].str.len_chars()
	deletion_mask = (mut_info["Ref"].str.len_chars() > 1) & (mut_info["Alt"].str.len_chars() == 1)
	insertion_mask = (mut_info["Ref"].str.len_chars() == 1) & (mut_info["Alt"].str.len_chars() > 1)
	
	mut_info = (
		mut_info.with_columns([
			pl.when(snv_mask)
			.then(mut_info["Ref"].str.len_chars().cast(pl.Utf8) + "-snv")
			.when(deletion_mask)
			.then((mut_info["Ref"].str.len_chars() - 1).cast(pl.Utf8) + "-del")
			.when(insertion_mask)
			.then((mut_info["Alt"].str.len_chars() - 1).cast(pl.Utf8) + "-ins")
			.otherwise(pl.lit("-"))
			.alias("Mut_type")
		])
	)


	# Add the simpleRepeat_TRF column by checking if the variant is in a simple repeat region
	# simple repeats bed obtained from ucsc table browser is used
	if simple_repeats_bed is not None:
		mut_info = annotate_simple_repeats(mut_info, simple_repeats_bed)
	else:
		# If no simple repeats bed provided, fill with "-"
		mut_info = mut_info.with_columns(
			pl.lit("-").alias("SimpleRepeat_TRF")
		)


	# Reoder columns accodign to MicroSEC documentation
	mut_info = mut_info.select([
		"Sample",
		"Mut_type",
		"Chr",
		"Pos",
		"Ref",
		"Alt",
		"SimpleRepeat_TRF",
		"Neighborhood_sequence"
	])

	return mut_info



## Main Processing

### setup paths and variables
sample_paths = pl.read_csv("annot/bam_vcf_path-with_variants.absolute.tsv", separator="\t")
simple_repeat = return_path_if_exists("ref/hg19/ucsc.hg19.simple-repeat.bed", abs=True)
genome = pysam.FastaFile("ref/hg19/ucsc.hg19.fasta")
outdir_root = "vcf-micr-svf"

# Read in the simple repeats bed file to annotate variants later
simple_repeats_bed = pl.read_csv(
	simple_repeat, 
	separator="\t", 
	has_header=False, 
	new_columns=["chrom", "start", "end", "type", "score"],
).select(["chrom", "start", "end"])



### Create mutation information files for MicroSEC
for i, path in enumerate(sample_paths["vcf_path"]):

	sample = sample_paths[i, "sample_name"]

	print(f"{i}. Creating MicroSEC input for {sample} from VCF: {path}")
	
	mut_info = make_mut_info(
		vcf_path=path, 
		genome=genome, 
		sample_name=sample, 
		simple_repeats_bed=simple_repeats_bed,
		ct_only=False, 
		snv_only=False
	)
 
	output_dir = f"{outdir_root}/{sample}/inputs"
	os.makedirs(output_dir, exist_ok=True)
	output_path = f"{output_dir}/{sample}.microsec.mut-info.tsv"
	mut_info.write_csv(output_path, separator="\t")



### Create Sample Info for MicroSEC (polars)
mut_info_paths = sorted(glob.glob(f"{outdir_root}/*/*/*.microsec.mut-info.tsv"))
mut_info_suffix = ".microsec.mut-info.tsv"

## Subset based on chunk of bam files available on disk
subset = sample_paths[1426:]

rows = []
for i, sample_name in enumerate(subset["sample_name"]):

	print(f"Processing sample {i+1}: {sample_name}")	

	mut_info_file = os.path.abspath(f"{outdir_root}/{sample_name}/inputs/{sample_name}{mut_info_suffix}")
	bam_path = subset[i, "bam_path"]
	ref_path = os.path.abspath("ref/hg19/ucsc.hg19.fasta")
	read_length = get_read_length(bam_path)
	genome_build = "hg19"
	## Using illumina adapters as foundation medicine uses the illumina platform
	adapter_1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
	adapter_2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

	rows.append({
		"sample_name": sample_name,
		"mutation_file": mut_info_file,
		"bam_file": bam_path,
		"read_length": read_length,
		"adapter_1": adapter_1,
		"adapter_2": adapter_2, # Can be None/NA
		"genome_build": genome_build, # e.g. "hg38" or "Human"
		"panel_name": "TOP",
		"ref_genome": ref_path, # Path to fasta if needed for CRAM
		"simple_repeat_bed": simple_repeat # Optional
	})

sample_info = pl.DataFrame(rows)	

# Write out individual sample info files
for i, sample_name in enumerate(sample_info.get_column("sample_name")):
	sample_info[i].write_csv(f"{outdir_root}/{sample_name}/inputs/{sample_name}.microsec.sample_info.tsv", separator="\t")


### Create MicroSEC execution scripts
sample_info_paths = [os.path.abspath(path) for path in glob.glob(f"{outdir_root}/*/*/*.microsec.sample_info.tsv")]
script_dir = f"{outdir_root}/script_microsec"
os.makedirs(script_dir, exist_ok=True)
microsec_script_path = return_path_if_exists(f"{outdir_root}/microsec.R", abs=True)

for path in sample_info_paths:

	sample_name = os.path.basename(path).replace(".microsec.sample_info.tsv", "")
	script_path = f"{script_dir}/{sample_name}.microsec.sh"
	microsec_outdir = os.path.abspath(f"{outdir_root}/{sample_name}")
	
	with open(script_path, "w") as file:
		file.writelines("#!/usr/bin/env bash\n")
		file.writelines(f"Rscript {microsec_script_path} --sample_info '{path}' --output_dir '{microsec_outdir}'\n")
	
