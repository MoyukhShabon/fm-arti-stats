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

def annotate_res(mobsnvf_res_path: str, annotated_snvs_path: str) -> pl.DataFrame:
	raw_res = pl.read_csv(mobsnvf_res_path, separator = "\t")
	if "is_vus" in raw_res.columns:
		return raw_res
	else:
		annotated_snvs = pl.read_csv(annotated_snvs_path, separator="\t")
		return raw_res.join(annotated_snvs, on = ["chrom", "pos", "ref", "alt"], how = "left")

# %% [markdown]
# ### VUS exclusion

no_bam_list = pl.read_csv("annot/vcf-no_bam.tsv", separator="\t")
xml_snv_annot_paths = sorted(glob.glob('xml-variants/*/*.with-indels.tsv'))


no_var_less_vus = []
no_vus = []

for i, path in enumerate(xml_snv_annot_paths):
	sample_name = path.split("/")[-2]
	
	if sample_name in no_bam_list["sample_name"].to_list():
		print(f"{i+1}. Sample {sample_name} has no associated BAM. skipping...")
		continue
	
	print(f"{i+1}. Excluding VUS from mobsnvf and microsec results of {sample_name}")
	
	xml_path = return_if_exists(f"data/{sample_name}.xml")
	ffpe_snvf_path = return_if_exists(f"xml-ffpe-snvf/{sample_name}/{sample_name}.mobsnvf.ffpe.snv")
	oxog_snvf_path = return_if_exists(f"xml-oxog-snvf/{sample_name}/{sample_name}.mobsnvf.oxog.snv")
	micr_svf_path = return_if_exists(f"xml-micr-svf/{sample_name}/{sample_name}.microsec.tsv")
	
	mobsnvf_ffpe_annot = annotate_res(ffpe_snvf_path, path)
	mobsnvf_oxog_annot = annotate_res(oxog_snvf_path, path)
	msec_micr_annot = annotate_res(micr_svf_path, path)

	mobsnvf_ffpe_no_vus = mobsnvf_ffpe_annot.filter(~pl.col("is_vus"))
	mobsnvf_oxog_no_vus = mobsnvf_oxog_annot.filter(~pl.col("is_vus"))
	msec_micr_no_vus = msec_micr_annot.filter(~pl.col("is_vus"))

	mobsnvf_ffpe_vus_only = mobsnvf_ffpe_annot.filter(pl.col("is_vus"))
	mobsnvf_oxog_vus_only = mobsnvf_oxog_annot.filter(pl.col("is_vus"))
	msec_micr_vus_only = msec_micr_annot.filter(pl.col("is_vus"))

	
	if mobsnvf_ffpe_no_vus.is_empty():
		print(f"\t{sample_name} has no variants after VUS exclusion")
		no_var_less_vus.append(sample_name)
	else:
		mobsnvf_ffpe_no_vus.write_csv(ffpe_snvf_path.replace(".snv", ".no_vus.snv"), separator="\t")
		mobsnvf_oxog_no_vus.write_csv(oxog_snvf_path.replace(".snv", ".no_vus.snv"), separator="\t")
		msec_micr_no_vus.write_csv(micr_svf_path.replace(".tsv", ".no_vus.tsv"), separator="\t")


	if mobsnvf_ffpe_vus_only.is_empty():
		print(f"\t{sample_name} has no VUS")
		no_vus.append(sample_name)
	else:
		mobsnvf_ffpe_vus_only.write_csv(ffpe_snvf_path.replace(".snv", ".vus_only.snv"), separator="\t")
		mobsnvf_oxog_vus_only.write_csv(oxog_snvf_path.replace(".snv", ".vus_only.snv"), separator="\t")
		msec_micr_vus_only.write_csv(micr_svf_path.replace(".tsv", ".vus_only.tsv"), separator="\t")

   
pl.DataFrame(no_var_less_vus).write_csv("annot/xml-no_snv-vus_filtration.txt", separator="\t", include_header=False)
pl.DataFrame(no_vus).write_csv("annot/xml-no_vus.txt", separator="\t", include_header=False)

