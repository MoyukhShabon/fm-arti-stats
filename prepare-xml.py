#!/usr/bin/env bash

import polars as pl
from lxml import etree
import os
import glob
import re

import sys

sys.path.append(os.path.abspath("common-ffpe-snvf/python"))

from mutation_signatures import read_variants

# %% [markdown]
# ## Functions

def str_to_bool(bool_str: str) -> bool | str:
	if bool_str.lower() == "true":
		return True
	elif bool_str.lower() == "false":
		return False
	else:
		return bool_str

def return_path_if_exists(path: str) -> str:
	if not os.path.exists(path):
		raise FileNotFoundError(f"File not found: {path}")
	else:
		return path


def snv_filter(df: pl.DataFrame) -> pl.DataFrame:
	# Filter for SNVs only: Length of Ref and Alt must be 1
	df = df.filter(
		(pl.col("ref").str.len_chars() == 1) & 
		(pl.col("alt").str.len_chars() == 1)
	)
	return df

def get_variant_properties(xml_root: etree._Element) -> pl.DataFrame:

	variant_properties = []

	for property in xml_root.xpath("//VariantProperty"):
		
		properties = {
			"is_vus" : str_to_bool(property.get("isVUS")),
			"gene": property.get("geneName"),
			"variant_name": property.get("variantName"),
		}
		
		variant_properties.append(properties)
		
	return pl.DataFrame(variant_properties)

def parse_variants_from_xml(xml_path: str, exclude_vus: bool = False, include_indels: bool = False) -> pl.DataFrame:
	"""
	Parses an XML report to extract variants and enriches them with REF/ALT alleles 
	from a companion VCF file located in the same directory.

	Args:
		xml_path: The file path to the input XML report.
		exclude_vus: If True, filters out variants marked as VUS.
		include_indels: If True, returns SNVs and INDELs. If False, filters for SNVs only.

	Returns:
		A Polars DataFrame containing the parsed variants.
	"""

	# ---- 1. Load and Prep VCF ----
	vcf_path = xml_path.replace(".xml", ".vcf")
	
	if not os.path.exists(vcf_path):
		raise FileNotFoundError(f"Companion VCF not found at {vcf_path}")

	# Read VCF, strip '#' from headers, lowercase columns
	vcf = pl.read_csv(
		vcf_path, 
		comment_prefix="##", 
		separator="\t", 
		infer_schema_length=1000,
		null_values=['.'],
		columns=["#CHROM", "POS", "REF", "ALT"]
	).rename(lambda x: x.strip("#").lower())

	# ---- 2. Parse XML Annotations ----
	try:
		tree = etree.parse(xml_path)
		root = tree.getroot()
	except etree.XMLSyntaxError as e:
		raise ValueError(f"Error parsing XML file: {e}")

	namespaces = {
		'vr': 'http://foundationmedicine.com/compbio/variant-report-external'
	}

	xml_variants_data = []

	# Use XPath to find all 'short-variant' elements
	for variant in root.xpath('.//vr:short-variant', namespaces=namespaces):
		cds_effect_str = variant.get('cds-effect')
		position_str = variant.get('position')
		equivocal_str = variant.get('equivocal')

		if not all([position_str, cds_effect_str]):
			continue

		try:
			# Parse chromosome and position from XML (e.g., "chr7:140453136")
			chrom, pos_str = position_str.split(':')
			pos = int(pos_str)

			variant_dict = {
				'chrom': chrom,
				'pos': pos,
				# 'ref' and 'alt' will come from VCF join
				'gene': variant.get('gene'),
				'depth': int(variant.get('depth')),
				'cds_effect': cds_effect_str,
				'protein_effect': variant.get('protein-effect'),
				'allele_fraction': float(variant.get('allele-fraction')),
				'functional_effect': variant.get('functional-effect'),
				'transcript': variant.get('transcript'),
				'strand': variant.get('strand'),
				'equivocal': str_to_bool(equivocal_str),
			}
			xml_variants_data.append(variant_dict)

		except (ValueError, IndexError) as e:
			print(f"\tSkipping variant due to parsing error: {etree.tostring(variant)}")
			continue

	if not xml_variants_data:
		print(f"\tNo variants detected in XML: {xml_path}")
		return pl.DataFrame(), pl.DataFrame()

	df_xml = pl.DataFrame(xml_variants_data)

	# ---- 3. Merge XML and VCF ----
	
	# Handle cases where VCF does not have variants present in XML. This should in theory never happen.
	var_missing_in_vcf = (
    	df_xml
    	.join(vcf, on=["chrom", "pos"], how="anti")
    	.with_columns(pl.lit(os.path.basename(xml_path).replace(".xml", "")).alias("sample_name"))
	)

	if var_missing_in_vcf.height > 0:
		print(f"\tWarning: {var_missing_in_vcf.height} variants from XML not found in VCF.")
	
	# Left join XML variants with VCF to get genomic REF/ALT
	# We join on Chrom and Pos. 
	df = df_xml.join(
		vcf, 
		on=["chrom", "pos"], 
		how="left"
	)

	# ---- 4. Indel Filtering ----
	
	if not include_indels:
		df = snv_filter(df)

	# ---- 5. VUS Annotation ----

	var_properties = get_variant_properties(root)
	
	if not var_properties.is_empty():
		var_properties = var_properties.filter(pl.col("is_vus")).select(["is_vus", "variant_name"])
		df = (
			df
			.join(var_properties, left_on="protein_effect", right_on="variant_name", how="left")
			.with_columns(pl.col("is_vus").fill_null(False))
		)
	else:
		df = df.with_columns(pl.lit(False).alias("is_vus"))

	# ---- 6. VUS Exclusion ----

	if exclude_vus:
		df = df.filter(~pl.col("is_vus"))
	
	# ---- 7. Formatting ----
	
	# Here we keep normalized numbers/letters but sort numerically.
	chrom_sort_key = (
		pl.when(pl.col("chrom") == "X").then(pl.lit(23, dtype=pl.Int64))
		.when(pl.col("chrom") == "Y").then(pl.lit(24, dtype=pl.Int64))
		.when((pl.col("chrom") == "M") | (pl.col("chrom") == "MT")).then(pl.lit(25, dtype=pl.Int64))
		.otherwise(pl.col("chrom").str.extract(r"(\d+)", 1).cast(pl.Int64))
		.fill_null(99) # Place any other contigs at the very end
	)
	
	schema = [
		'chrom', 'pos', 'ref', 'alt', 
		'gene', 'is_vus', 'depth', 'cds_effect', 'protein_effect',
		'allele_fraction', 'functional_effect', 
		'transcript', 'strand', 'equivocal',
	]
	
	# Select columns, ensuring they exist
	final_cols = [c for c in schema if c in df.columns]
	df = df.select(final_cols).sort(chrom_sort_key, pl.col("pos"))
	
	return df, var_missing_in_vcf

# %% [markdown]
# ## Parse XML to SNV

outdir_root = "xml-snvs"

xml_paths = sorted(glob.glob("data/*.xml"))

xml_table = pl.DataFrame({
	"sample_name": [os.path.basename(p).replace(".xml", "") for p in xml_paths],
	"xml_path": [os.path.abspath(p) for p in xml_paths]
})

bam_vcf_table = pl.read_csv("annot/bam_vcf_path.absolute.tsv", separator="\t")

bam_xml_table = bam_vcf_table.join(xml_table, on="sample_name", how="inner")
bam_xml_table.write_csv("annot/bam_xml_path.absolute.tsv", separator="\t")

no_variants = []
missing_vars = []

for i, path in enumerate(xml_paths):
	
	print(f"{i+1}. Processing {path}")	
 
	sample_name = os.path.basename(path).replace(".xml", "")
	outpath = f"{outdir_root}/{sample_name}"
	os.makedirs(outpath, exist_ok=True)
 
	variants, missing_in_vcf = parse_variants_from_xml(path, include_indels=True)
	missing_vars.append(missing_in_vcf)
	
	if variants.is_empty():
		print(f'Dataframe for "{sample_name}" is empty. Skipping writing to disk')
		no_variants.append(
			pl.DataFrame({"sample_name": sample_name, "xml_path": path})
		)
		continue
 
	variants.write_csv(f"{outpath}/{sample_name}.with-indels.tsv", separator="\t")
	snv_filter(variants).select(["chrom", "pos", "ref", "alt"]).write_csv(f"{outpath}/{sample_name}.snv", separator="\t")

no_variants = pl.concat(no_variants)
missing_vars = pl.concat(missing_vars, how="diagonal_relaxed")
no_variants.write_csv("annot/xml-no_snv.tsv", include_header=False, separator="\t")
missing_vars.write_csv("annot/xml-variants_missing_in_vcf.tsv", separator="\t")

with_snv_with_bam = bam_xml_table.join(no_variants, on="sample_name", how="anti")
with_snv_with_bam.write_csv("annot/bam_xml_path-with_variants.absolute.tsv", separator="\t")

xml_no_snv_no_bam = xml_table.join(with_snv_with_bam, on="sample_name", how="anti")
xml_no_snv_no_bam.write_csv("annot/xml-no_snv-no_bam.tsv", include_header=False, separator="\t")

# %% [markdown]
# ## Subset XML variants
# Make a subset of the VCF SNVF results to just retain the variants in the XML. Variant in the XML are a subset of the variants in the VCF

no_bam_list = pl.read_csv("annot/vcf-no_bam.tsv", separator="\t")["sample_name"].to_list()

xml_snv_paths = sorted(glob.glob("xml-snvs/*/*.with-indels.tsv"))
mobsnvf_paths = sorted(glob.glob("vcf-ffpe-snvf/*/*.mobsnvf.ffpe.snv") + glob.glob("vcf-oxog-snvf/*/*.mobsnvf.oxog.snv"))

def subset_variants(mobsnvf_res_path: str, xml_snvs: pl.DataFrame, write_output: bool = False, annotations: bool = False) -> pl.DataFrame:
	mobsnvf_res = read_variants(mobsnvf_res_path).rename({"fobp":"FOBP"})

	if annotations:
		mobsnvf_res_subset =  mobsnvf_res.join(xml_snvs, on = ["chrom", "pos", "ref", "alt"], how="inner")
	else:
		mobsnvf_res_subset =  mobsnvf_res.join(xml_snvs, on = ["chrom", "pos", "ref", "alt"], how="semi")

	if write_output:
		outpath = mobsnvf_res_path.replace("vcf", "xml")
		os.makedirs(os.path.dirname(outpath), exist_ok=True)
		mobsnvf_res_subset.write_csv(outpath, separator="\t")

		print(f"\tXML subset written to : {outpath}")

	return mobsnvf_res_subset

for i, path in enumerate(xml_snv_paths):
	
	sample_name = os.path.basename(path).split(".")[0]
	print(f"{i+1}. Subsetting VCF SNVF results with variants in: {path}")

	sample_mobsnvf_path = [path for path in mobsnvf_paths if sample_name in path]

	if not sample_mobsnvf_path:
		message = f"\tVCF SNVF results missing for {sample_name}"
		if sample_name in no_bam_list:
			message = f"{message} as it has no associated BAM"

		print(message)
		continue
	
	xml_snvs = read_variants(path, sort_variants=True)

	for mobsnvf_res_path in sample_mobsnvf_path:
		subset_variants(mobsnvf_res_path, xml_snvs, write_output=True, annotations=True)

