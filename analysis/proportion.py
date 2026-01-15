#!/usr/bin/env python
import polars as pl
import os
import glob
from lxml import etree

# %% [markdown]
# ## Functions

def get_test_type(xml_path: int) -> str:
	
	if not os.path.exists(xml_path):
		raise FileNotFoundError(f"No file found at {xml_path}")
	
	tree = etree.parse(xml_path)
	root = tree.getroot()
	return (
		root
		.xpath("//TestType")[0]
		.text
		.replace(" ", "")
		.replace("FoundationOneCDX", "FoundationOneCDx")
		.replace("FoundationOneLiquidCDx", "FoundationOneLiquidDx")
	)


def calc_arti_prop(artifact_pred_paths: list, arti_type: str, variant_source: str) -> pl.DataFrame:
	arti_prop = []

	# c = 1
	for path in artifact_pred_paths:
		
		artifact_pred = pl.read_csv(path, separator="\t", infer_schema_length=1000)

		if artifact_pred.is_empty():
			# print(f"{c}. Skipping {path} as it is empty. \n\tLikely Reasons include no corresponding BAMs, no variants in the VCF, or no variants after VUS exclusion")
			# c += 1
			continue

		variants = artifact_pred.shape[0]
		
		if arti_type.lower() == "micr":
			artifacts = artifact_pred.filter(~pl.col("msec_filter_1234").is_null()).shape[0]
		else:
			artifacts = artifact_pred.filter(~pl.col("pred")).shape[0]

		sample_name = os.path.basename(path).split(".")[0]
		xml_path = f"../data/{sample_name}.xml"
		test_type = get_test_type(xml_path)

		sample_stats = {
			"sample_name" : sample_name,
			"source": variant_source.upper(),
			"test_type": test_type,
			"arti_type": arti_type.upper(),
			"n_variants" : variants,
			"n_real" : (variants - artifacts),
			"n_artifacts" : artifacts,
			"proportion" : (artifacts / variants)
		}

		arti_prop.append(sample_stats)

	arti_prop = pl.DataFrame(arti_prop).sort("proportion", descending=True)
	return arti_prop

def get_stats(arti_prop: pl.DataFrame, damage_type: str = None, variant_source: str = None) -> None:
	
	prop_test_type = (
		arti_prop
		.group_by("test_type")
		.agg(
			pl.col("proportion").mean().alias("mean_proportion"),
			pl.col("proportion").count().alias("count"),
			pl.col("proportion").max().alias("max_proportion")
		)
	)

	prop_test_type_with_artifact = (
		arti_prop
		.filter(pl.col("proportion") > 0)
		.group_by("test_type")
		.agg(
			pl.col("proportion").mean().alias("mean_proportion"),
			pl.col("proportion").count().alias("count_>=1_artifact")
		)
	)
	
	mean_pct = arti_prop["proportion"].mean() * 100
	max_pct = arti_prop["proportion"].max() * 100
	n_samples = arti_prop.shape[0]
	has_artifact_df = arti_prop.filter(pl.col("proportion") > 0)
	n_with_artifacts = has_artifact_df.shape[0]
	mean_pct_arti_samples = has_artifact_df["proportion"].mean() * 100

	arti_type = damage_type.upper()
	print("Damage Type:", arti_type)
	print("Variant Source:", variant_source.upper())

	print(f"Samples analyzed: {n_samples}")
	print(f"Samples with ≥1 predicted {arti_type} artifact: {n_with_artifacts} ({n_with_artifacts / n_samples * 100:.1f}%)")
	print(f"Mean proportion of {arti_type} artifactual SNVs: {mean_pct:.2f}%")
	print(f"Max proportion of {arti_type} artifactual SNVs: {max_pct:.2f}%")
	print(f"Mean proportion of {arti_type} artifactual SNVs within samples with ≥1 predicted artifact: {mean_pct_arti_samples:.2f}%")
 
	print("\nResults stratified by test type:\n")
	for i in range(prop_test_type.shape[0]):

		test_type = prop_test_type[i, "test_type"]

		print(f"\t{test_type} samples analyzed: {prop_test_type[i, "count"]}")
		print(f"\t{test_type} max {arti_type} artifact proportion within samples: {prop_test_type[i, "max_proportion"] * 100:.2f}%")
		print(f"\t{test_type} mean {arti_type} artifact proportion within samples: {prop_test_type[i, "mean_proportion"] * 100:.2f}%")

		prop_test_type_with_artifact_filtered = prop_test_type_with_artifact.filter(pl.col("test_type") == test_type)

		if prop_test_type_with_artifact_filtered.shape[0]:
			print(f"\t{test_type} samples with {arti_type} artifact: {prop_test_type_with_artifact_filtered[0, "count_>=1_artifact"]}")
			print(f"\t{test_type} mean artifact proportion within samples with ≥1 detected {arti_type} artifacts: {prop_test_type_with_artifact_filtered[0, "mean_proportion"] * 100:.2f}%")

		print()
	print("-----------------\n")

def get_res(variant_source: str, damage_type: str, fp_cut: float, write_data: bool = True, exclude_vus = False) -> pl.DataFrame:
	
	if not exclude_vus:
		vus_token = ""
	else:
		vus_token = ".no_vus"
		print("VUS Excluded")

	if damage_type.lower() == "micr":
		search_pattern = f"../{variant_source.lower()}-{damage_type.lower()}-svf/*/*.microsec{vus_token}.tsv"
	else:
		search_pattern = f"../{variant_source.lower()}-{damage_type.lower()}-snvf/*/*.mobsnvf.{damage_type.lower()}{vus_token}.pred_fp-cut_{fp_cut:.0e}.tsv"
	
	arti_pred_paths = glob.glob(search_pattern)
	
	if not arti_pred_paths:
		raise FileNotFoundError (f"No predictions found at {search_pattern}")

	proportions = calc_arti_prop(arti_pred_paths, damage_type, variant_source)

	if fp_cut:
		print(f"FP-cut: {fp_cut}")

	get_stats(proportions, damage_type=damage_type, variant_source=variant_source)

	if write_data:
		if damage_type.lower() == "micr":
			proportions.write_csv(f"{damage_type.lower()}_proportions_per_sample.{variant_source.lower()}.tsv", separator="\t")
		else:
			proportions.write_csv(f"{damage_type.lower()}_proportions_per_sample.{variant_source.lower()}{vus_token}.fp-cut_{fp_cut:.0e}.tsv", separator="\t")
	
	return proportions

# %% [markdown]
# ## MICR Artifacts

get_res("VCF", "MICR", None)

get_res("XML", "MICR", None)

get_res("XML", "MICR", None, exclude_vus=True)

# %% [markdown]
# ## FP-CUT=1e-08

# %% [markdown]
# ### FFPE Artifacts

get_res("VCF", "FFPE", 1e-08)

get_res("XML", "FFPE", 1e-08)

get_res("XML", "FFPE", 1e-08, exclude_vus=True)

# %% [markdown]
# ### OXOG Artifacts

get_res("VCF", "OXOG", 1e-08)

get_res("XML", "OXOG", 1e-08)

get_res("XML", "OXOG", 1e-08, exclude_vus=True)

# %% [markdown]
# ## FP-CUT = 5e-01
# 
# This is not a good cut off according to our testing.

# %% [markdown]
# ### FFPE artifacts

get_res("vcf", "ffpe", 5e-01)

get_res("xml", "ffpe", 5e-01)

get_res("xml", "ffpe", 5e-01, exclude_vus=True)

# %% [markdown]
# ### OXOG Artifacts

get_res("vcf", "oxog", 5e-01)

get_res("xml", "oxog", 5e-01)

get_res("xml", "oxog", 5e-01, exclude_vus=True)


