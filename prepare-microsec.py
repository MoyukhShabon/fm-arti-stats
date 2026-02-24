#!/usr/bin/env python
import polars as pl
import os
import glob
import pysam
import sys

sys.path.append("common-ffpe-snvf/python")
from microsec_preparation import *
from common import return_path_if_exists


sample_paths = pl.read_csv("annot/bam_vcf_path-with_variants.absolute.tsv", separator="\t")
simple_repeat = return_path_if_exists("ref/hg19/ucsc.hg19.simple-repeat.bed", abs=True)
ref_path = "ref/hg19/ucsc.hg19.fasta"
outdir_root = "vcf-micr-svf"

# Read in the simple repeats bed file to annotate variants later
simple_repeats_bed = pl.read_csv(
	simple_repeat, 
	separator="\t", 
	has_header=False, 
	new_columns=["chrom", "start", "end", "type", "score"],
).select(["chrom", "start", "end"])


### Create Mutation Info for MicroSEC (polars)
for i, vcf_path in enumerate(sample_paths["vcf_path"]):

	sample = sample_paths[i, "sample_name"]

	print(f"{i}. Creating MicroSEC input for {sample} from VCF: {vcf_path}")
	
	mut_info = make_mut_info(
		vcf_path=vcf_path, 
		ref_fasta_path=ref_path, 
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
subset = sample_paths[:]

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
	
