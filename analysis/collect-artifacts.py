#!/usr/bin/env python
import glob
import os
import polars as pl
from lxml import etree

## TODO: Implement functionality to collect all variants regardless of artifact or not

# %% [markdown]
# ### Useful Setup

# Useful column names and order
schema_mobsnvf = [
	'sample_name',
	'chrom',
	'pos',
	'ref',
	'alt',
	'arti_type',
	'source',
	'test_type',
	'is_vus',
	'gene',
	'depth',
	'cds_effect',
	'protein_effect',
	'allele_fraction',
	'functional_effect',
	'transcript',
	'strand',
	'equivocal',
	'score',
	'q',
	'pred',
]

schema_msec = [
    'sample_name', 
    # 'mut_type', 
    'chrom', 
    'pos', 
    'ref', 
    'alt', 
    'arti_type',
    'source',
	'test_type',
	'is_vus',
	'gene',
	'depth',
	'cds_effect',
	'protein_effect',
	'allele_fraction',
	'functional_effect',
	'transcript',
	'strand',
	'equivocal',
    # 'simplerepeat_trf', 
    # 'neighborhood_sequence', 
    # 'read_length', 
    # 'total_read', 
    # 'soft_clipped_read', 
    # 'flag_hairpin', 
    # 'pre_support_length', 
    # 'post_support_length', 
    # 'short_support_length', 
    # 'pre_farthest', 
    # 'post_farthest', 
    # 'low_quality_base_rate_under_q18', 
    # 'low_quality_pre', 
    # 'low_quality_post', 
    # 'distant_homology_rate', 
    # 'soft_clipped_rate', 
    # 'prob_filter_1', 
    # 'prob_filter_3_pre', 
    # 'prob_filter_3_post', 
    'filter_1_mutation_intra_hairpin_loop', 
    'filter_2_hairpin_structure', 
    'filter_3_microhomology_induced_mutation', 
    'filter_4_highly_homologous_region', 
    'filter_5_soft_clipped_reads', 
    'filter_6_simple_repeat', 
    'filter_7_mutation_at_homopolymer', 
    'filter_8_low_quality', 
    'msec_filter_123', 
    'msec_filter_1234', 
    'msec_filter_all', 
    # 'comment', 
]

# Chromosome string to numeric key. Used for sorting
chrom_indices = {
    'chr1': 1, 'chr2': 2, 'chr3': 3, 'chr4': 4, 'chr5': 5,
    'chr6': 6, 'chr7': 7, 'chr8': 8, 'chr9': 9, 'chr10': 10,
    'chr11': 11, 'chr12': 12, 'chr13': 13, 'chr14': 14, 'chr15': 15,
    'chr16': 16, 'chr17': 17, 'chr18': 18, 'chr19': 19, 'chr20': 20,
    'chr21': 21, 'chr22': 22, 'chrX': 23, 'chrY': 24,
    'chrM': 25, 'chrMT': 25  # Both M and MT map to the same index
}

# %% [markdown]
# ### Helper Functions

def annotate_variants(annot_path: str, variant_df:pl.DataFrame) -> pl.DataFrame:
	annot = pl.read_csv(annot_path, separator="\t")
	return variant_df.join(annot, on = ["chrom", "pos", "ref", "alt"], how="left")

def get_test_type(xml_path: int) -> str:
	
	if not os.path.exists(xml_path):
		raise FileNotFoundError(f"No file found at {xml_path}")
	
	tree = etree.parse(xml_path)
	root = tree.getroot()
	return root.xpath("//TestType")[0].text

def snv_filter(df):
	return (
		df.filter(
			(pl.col("ref").str.len_chars() == 1),
			(pl.col("alt").str.len_chars() == 1)
		)
	)

# %% [markdown]
# ### MOBSNVF | FFPE and OXOG Artifact Collection

fp_cut = 1e-08

# Path to MOBSNVF predictions
pred_paths = sorted(
    glob.glob(f"../*/*/*.mobsnvf.ffpe.pred_fp-cut_{fp_cut:.0e}.tsv") +
	glob.glob(f"../*/*/*.mobsnvf.oxog.pred_fp-cut_{fp_cut:.0e}.tsv")
)

collected_artifacts = []

for i, path in enumerate(pred_paths):
	
	sample_name = os.path.basename(path).split(".")[0]
	arti_type = os.path.basename(path).split(".")[2].upper()
	source = path.split("/")[1].split("-")[0].upper()
	test_type = get_test_type(f"../data/{sample_name}.xml")
	
	print(f"{i+1}. Processing {path}")
	
	pred = (
		pl.read_csv(path, separator="\t")
		.with_columns(
			pl.lit(sample_name).alias("sample_name"), 
			pl.lit(arti_type).alias("arti_type"),  
			pl.lit(source).alias("source")
		)
	)

	arti = (
		pred
		.filter(~pl.col("pred"))
		.with_columns(pl.lit(test_type).alias("test_type"))
	)
 
	collected_artifacts.append(arti)
	

ffpe_oxog_arti = (
    pl.concat(collected_artifacts, how="diagonal_relaxed")
	# Remove duplicated artifacts i.e artifacts that appear both in XML and VCF. Only the XML artifacts are retained
	.filter(~((pl.struct(["sample_name", "chrom", "pos", "ref", "alt"]).is_duplicated()) & (pl.col("source") == "VCF")))
	# Reorder and select necessary columns
	.select(schema_mobsnvf)
    .rename({
        "score" : "mobsnvf_score",
        "q": "mobsnvf_q",
        "pred": "mobsnvf_pred"
	})
)

# ffpe_oxog.write_csv("all_artifacts.mobsnvf.pred_fp-cut_1e-08.tsv", separator="\t")

# %% [markdown]
# ## MicroSEC | MICR artifact collection

microsec_paths = sorted(
    glob.glob("../*/*/*.microsec.tsv")
)

msec_collected = []

for i, path in enumerate(microsec_paths):
	sample_name = os.path.basename(path).split(".")[0]
	arti_type = "MICR"
	source = path.split("/")[1].split("-")[0].upper()
	test_type = get_test_type(f"../data/{sample_name}.xml")
	
	
	print(f"{i+1}. Processing {path}")
	
	msec_sample = (
		pl.read_csv(path, separator="\t", infer_schema_length=1000)
		.filter(pl.col("msec_filter_1234").is_not_null())
		.with_columns(
			pl.lit(arti_type).alias("arti_type"),  
			pl.lit(source).alias("source")
		)
		.rename(lambda x : x.lower())
		.rename({"sample": "sample_name"})
		.with_columns(pl.lit(test_type).alias("test_type"))
	)

	if "chr" in msec_sample.columns:
		msec_sample = msec_sample.rename({"chr" : "chrom"})

	msec_collected.append(msec_sample)

micr_arti = (
	pl.concat(msec_collected, how="diagonal_relaxed")
	# Reorder and select necessary columns
	.select(schema_msec)
	# Only the XML artifacts are retained if they also appear in the VCF
    .filter(~((pl.struct(["sample_name", "chrom", "pos", "ref", "alt"]).is_duplicated()) & (pl.col("source") == "VCF")))
    .rename(lambda x: x.replace("msec_", "").replace("filter_", "msec_filter_"))
)


# %% [markdown]
# ### Combine Results

all_artifacts = (
	pl.concat([ffpe_oxog_arti, micr_arti], how="diagonal_relaxed")
	# Sort
	.with_columns(pl.col("chrom").map_elements(lambda x : chrom_indices.get(x, 99), return_dtype=int).alias("chrom_n"))
	.sort(["sample_name", "chrom_n", "pos", "source", "arti_type"])
	.drop("chrom_n")
)

all_artifacts.write_csv("all-artifacts_microsec_mobsnvf_fp-cut_1e-08.tsv", separator="\t")
all_artifacts.pipe(snv_filter).write_csv("all-artifacts_microsec_mobsnvf_fp-cut_1e-08.snv-only.tsv", separator="\t")

