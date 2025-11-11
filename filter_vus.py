#!/usr/bin/env bash

import polars as pl
from lxml import etree
import os
import glob
import re

# %% [markdown]
# ### Functions

def str_to_bool(bool_str: str) -> bool | str:
	if bool_str.lower() == "true":
		return True
	elif bool_str.lower() == "false":
		return False
	else:
		return bool_str

def return_if_exists(path: str) -> str:
	if not os.path.exists(path):
		raise FileNotFoundError(f"File not found: {path}")
	else:
		return path

def get_variant_properties(xml_path: str) -> pl.DataFrame:

	tree = etree.parse(xml_path)
	root = tree.getroot()

	variant_properties = []

	for property in root.xpath("//VariantProperty"):
		
		properties = {
			"is_vus" : str_to_bool(property.get("isVUS")),
			"gene": property.get("geneName"),
			"variant_name": property.get("variantName"),
		}
		
		variant_properties.append(properties)
		
	return pl.DataFrame(variant_properties)

def annotate_mobsnvf_res(mobsnvf_res_path: str, annotated_snvs: pl.DataFrame) -> pl.DataFrame:
	mobsnvf_raw = pl.read_csv(mobsnvf_res_path, separator = "\t")
	if len(mobsnvf_raw.columns) > 5:
		return mobsnvf_raw
	else:
		return mobsnvf_raw.join(annotated_snvs, on = ["chrom", "pos", "ref", "alt"], how = "left")

# %% [markdown]
# ### VUS exclusion

no_bam_list = pl.read_csv("annot/vcf-no_bam.tsv", separator="\t")
xml_snv_annot_paths = sorted(glob.glob('xml-snvs/*/*.tsv'))


no_variants = []

for i, path in enumerate(xml_snv_annot_paths):
	sample_name = path.split("/")[-2]
	
	if sample_name in no_bam_list["sample_name"].to_list():
		print(f"{i+1}. Sample {sample_name} has no associated BAM. skipping...")
		continue
	
	print(f"{i+1}. Excluding VUS from mobsnvf results of {sample_name}")
	
	xml_path = return_if_exists(f"data/{sample_name}.xml")
	ffpe_snvf_path = return_if_exists(f"xml-ffpe-snvf/{sample_name}/{sample_name}.mobsnvf.ffpe.snv")
	oxog_snvf_path = return_if_exists(f"xml-oxog-snvf/{sample_name}/{sample_name}.mobsnvf.oxog.snv")
	
	annotated_snvs = pl.read_csv(path, separator="\t")
	mobsnvf_ffpe_annot = annotate_mobsnvf_res(ffpe_snvf_path, annotated_snvs)
	mobsnvf_oxog_annot = annotate_mobsnvf_res(ffpe_snvf_path, annotated_snvs)
	
	if ("is_vus" not in mobsnvf_ffpe_annot.columns):
		# XMLs seem to only have variant properties with VUS being true
		# However a VUS == TRUE filter is applied just in case
		var_properties_vus = get_variant_properties(xml_path).filter(pl.col("is_vus"))
		
		# Subset variants
		mobsnvf_ffpe_no_vus = mobsnvf_ffpe_annot.join(var_properties_vus, left_on="protein_effect", right_on="variant_name", how="anti")
		mobsnvf_oxog_no_vus = mobsnvf_oxog_annot.join(var_properties_vus, left_on="protein_effect", right_on="variant_name", how="anti")

	else:
		mobsnvf_ffpe_no_vus = mobsnvf_ffpe_annot.filter(~pl.col("is_vus"))
		mobsnvf_oxog_no_vus = mobsnvf_oxog_annot.filter(~pl.col("is_vus"))
	
	
	keep_annotations = True
	
	if not keep_annotations:
		mobsnvf_ffpe_no_vus = mobsnvf_ffpe_no_vus.select("chrom", "pos", "ref", "alt", "FOBP")
		mobsnvf_oxog_no_vus = mobsnvf_oxog_no_vus.select("chrom", "pos", "ref", "alt", "FOBP")
	
	if mobsnvf_ffpe_no_vus.is_empty():
		print(f"\t{sample_name} has no variants after VUS exclusion")
		no_variants.append(sample_name)
		continue
	
	mobsnvf_ffpe_no_vus.write_csv(ffpe_snvf_path.replace(".snv", ".no_vus.snv"), separator="\t")
	mobsnvf_oxog_no_vus.write_csv(oxog_snvf_path.replace(".snv", ".no_vus.snv"), separator="\t")

   
pl.DataFrame(no_variants).write_csv("annot/xml-no_snv-vus_filtration.txt", separator="\t", include_header=False)
	


