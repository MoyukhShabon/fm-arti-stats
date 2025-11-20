#!/usr/bin/env python
import glob
import os
import polars as pl
from lxml import etree
import argparse

# %% [markdown]
# ### Parsed Arguments

parser = argparse.ArgumentParser(description="Collect artifacts that was predicted upstream with an specified fp-cut")
parser.add_argument("--fp-cut", "-c", default=1e-08, help="The False Positive Cutoff (fp-cut) value which was used to make predictions")
args = parser.parse_args()

# %% [markdown]
# ### Useful Stuff

# Useful column names and order
schema = [
	'sample_name',
	'source',
	'chrom',
	'pos',
	'ref',
	'alt',
	'arti_type',
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
# ### Functions

def annotate_variants(annot_path: str, variant_df:pl.DataFrame) -> pl.DataFrame:
	annot = pl.read_csv(annot_path, separator="\t")
	return variant_df.join(annot, on = ["chrom", "pos", "ref", "alt"], how="left")

def get_test_type(xml_path: int) -> str:
	
	if not os.path.exists(xml_path):
		raise FileNotFoundError(f"No file found at {xml_path}")
	
	tree = etree.parse(xml_path)
	root = tree.getroot()
	return root.xpath("//TestType")[0].text


# %% [markdown]
# ### Artifact Collection

fp_cut = float(args.fp_cut)

print(f"Compiling artifacts predicted with fp.cut = {fp_cut:.0e}\n")

# Path to MOBSNVF predictions
pred_paths = glob.glob(f"../*/*/*.mobsnvf.ffpe.pred_fp-cut_{fp_cut:.0e}.tsv") + glob.glob(f"../*/*/*.mobsnvf.oxog.pred_fp-cut_{fp_cut:.0e}.tsv")


all_artifacts = []

for i, path in enumerate(pred_paths):
	
	sample_name = os.path.basename(path).split(".")[0]
	arti_type = os.path.basename(path).split(".")[2].upper()
	source = path.split("/")[1].split("-")[0].upper()
	
	print(f"{i+1}. Processing {path}")
	
	pred = pl.read_csv(path, separator="\t").with_columns(pl.lit(sample_name).alias("sample_name"), pl.lit(arti_type).alias("arti_type"),  pl.lit(source).alias("source"))

	# # Annotate VCF SNVs where possible using XML
	# if source == "VCF":
	# 	annot_path = f"../xml-snvs/{sample_name}/{sample_name}.tsv"
	# 	if os.path.exists(annot_path):
	# 		pred = annotate_variants(annot_path, pred)

	arti = pred.filter(~pl.col("pred"))
 
	all_artifacts.append(arti)
	
all_artifacts = pl.concat(all_artifacts, how="diagonal_relaxed")

# Remove duplicated artifacts i.e artifacts that appear both in XML and VCF. Only the XML artifacts are retained
all_artifacts = all_artifacts.filter(~((all_artifacts[["sample_name", "chrom", "pos", "ref", "alt"]].is_duplicated()) & (all_artifacts["source"] == "VCF"))).sort("sample_name")

# Annotate Test Type
all_artifacts = all_artifacts.with_columns(pl.col("sample_name").map_elements(lambda x: get_test_type(f"../data/{x}.xml"), return_dtype=str).alias("test_type"))

# Reorder and select necessary columns
all_artifacts = all_artifacts.select(schema)

# Sort
all_artifacts = (
	all_artifacts
	.with_columns(pl.col("chrom").map_elements(lambda x : chrom_indices.get(x, 99), return_dtype=int).alias("chrom_n"))
	.sort(["sample_name", "chrom_n", "pos", "source", "arti_type"])
	.drop("chrom_n")
)

outdir = f"all_artifacts.mobsnvf.pred_fp-cut_{fp_cut:.0e}.tsv"
all_artifacts.write_csv(outdir, separator="\t")

print(f"\nTable with suspected artifacts written to: {outdir}")
