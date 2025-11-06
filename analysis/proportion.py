#!/usr/bin/env python
import polars as pl
import os
import glob
from lxml import etree

# %% [markdown]
# ## Functions

# %%
def get_test_type(xml_path: int) -> str:
	
	if not os.path.exists(xml_path):
		raise FileNotFoundError(f"No file found at {xml_path}")
	
	tree = etree.parse(xml_path)
	root = tree.getroot()
	return root.xpath("//TestType")[0].text


# %%
def calc_arti_prop(artifact_pred_paths: list) -> pl.DataFrame:
	arti_prop = []

	# c = 1
	for path in artifact_pred_paths:
		
		artifact_pred = pl.read_csv(path, separator="\t")

		if artifact_pred.is_empty():
			# print(f"{c}. Skipping {path} as it is empty. \n\tLikely Reasons include no corresponding BAMs, no variants in the VCF, or no variants after VUS exclusion")
			# c += 1
			continue

		snvs = artifact_pred.shape[0]
		
		artifacts = artifact_pred.filter(~pl.col("pred")).shape[0]

		sample_name = os.path.basename(path).split(".")[0]
		xml_path = f"../data/{sample_name}.xml"
		test_type = get_test_type(xml_path)

		sample_stats = {
			"sample_name" : sample_name,
			"test_type": test_type,
			"n_snv" : snvs,
			"n_real" : (snvs - artifacts),
			"n_artifacts" : artifacts,
			"proportion" : (artifacts / snvs)
		}

		arti_prop.append(sample_stats)

	arti_prop = pl.DataFrame(arti_prop).sort("proportion", descending=True)
	return arti_prop

# %%
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

	if damage_type:
		print("Damage Type:", damage_type.upper())

	if variant_source:
		print("Variant Source:", variant_source.upper())

	print(f"Samples analyzed: {n_samples}")
	print(f"Samples with ≥1 predicted artifact: {n_with_artifacts} ({n_with_artifacts / n_samples * 100:.1f}%)")
	print(f"Mean proportion of artifactual SNVs: {mean_pct:.2f}%")
	print(f"Max proportion of artifactual SNVs: {max_pct:.2f}%")
	print(f"Mean proportion of artifactual SNVs within samples with ≥1 predicted artifact: {mean_pct_arti_samples:.2f}%")
 
	print("\nResults stratified by test type:\n")
	for i in range(prop_test_type.shape[0]):

		test_type = prop_test_type[i, "test_type"]

		print(f"\t{test_type} samples: {prop_test_type[i, "count"]}")
		print(f"\t{test_type} max artifact proportion within samples: {prop_test_type[i, "max_proportion"] * 100:.2f}%")
		print(f"\t{test_type} mean artifact proportion within samples: {prop_test_type[i, "mean_proportion"] * 100:.2f}%")

		prop_test_type_with_artifact_filtered = prop_test_type_with_artifact.filter(pl.col("test_type") == test_type)

		if prop_test_type_with_artifact_filtered.shape[0]:
			print(f"\t{test_type} samples with >1 artifact: {prop_test_type_with_artifact_filtered[0, "count_>=1_artifact"]}")
			print(f"\t{test_type} mean artifact proportion within samples with ≥1 detected artifacts: {prop_test_type_with_artifact_filtered[0, "mean_proportion"] * 100:.2f}%")

		print()
	print("-----------------\n")

# %%
def get_res(variant_source: str, damage_type: str, fp_cut: float, write_data: bool = True, exclude_vus = False) -> pl.DataFrame:
	
	if not exclude_vus:
		vus_token = ""
	else:
		vus_token = ".no_vus"

	search_pattern = f"../{variant_source.lower()}-{damage_type.lower()}-snvf/*/*.mobsnvf.{damage_type.lower()}{vus_token}.pred_fp-cut_{fp_cut:.0e}.tsv"
	arti_pred_paths = glob.glob(search_pattern)
 
	if not arti_pred_paths:
		raise FileNotFoundError (f"No predictions made with fp-cut: {fp_cut:.0e} at {search_pattern}")

	proportions = calc_arti_prop(arti_pred_paths)

	print(f"FP-cut: {fp_cut}")
	get_stats(proportions, damage_type=damage_type, variant_source=variant_source)

	if write_data:
		proportions.write_csv(f"{damage_type}_proportions_per_sample.{variant_source.lower()}{vus_token}.fp-cut_{fp_cut:.0e}.tsv", separator="\t")
	
	return proportions

# %% [markdown]
# ## FFPE-VCF

# %%
get_res("vcf", "ffpe", 1e-08)

# %%
get_res("vcf", "ffpe", 5e-01)

# %% [markdown]
# #### Observations - FFPE
# 
# **2025/10/24**
# 
# Out of 204 samples analyzed so far.
# 
# - Using a stringent fp-cut of 1e-08: 
# 	- 39/204 samples (19.1%) was observed to possess artifacts. 
# 	- In terms of proportion of artifactual SNVs, the max was 25%, mean was 0.94 %
# 
# - Using a relaxed fp-cut of 5e-01: 
# 	- 9/204 samples (4.4%) was observed to possess artifacts. 
# 	- In terms of proportion of artifactual SNVs, the max was 25%, mean was 0.24 %
# 
# 
# **2025/10/27**
# 
# Out of 425 samples analyzed so far.
# 
# - Using a stringent fp-cut of 1e-08: 
# 	- 81/425 samples (19.1%) was observed to possess artifacts. 
# 	- In terms of proportion of artifactual SNVs, the max was 25%, mean was 0.93 %
# 
# - Using a relaxed fp-cut of 5e-01: 
# 	- 24/425 samples (5.6%) was observed to possess artifacts. 
# 	- In terms of proportion of artifactual SNVs, the max was 25%, mean was 0.26 %

# %% [markdown]
# ## OxoG-VCF

# %%
get_res("vcf", "oxog", 1e-08)

# %%
get_res("vcf", "oxog", 5e-01)

# %% [markdown]
# #### Observations - OxoG
# 
# **2025/10/24**
# 
# Out of 204 samples analyzed so far.
# 
# - Using a stringent fp-cut of 1e-08: 
# 	- 15/204 samples (7.3%) was observed to possess artifacts. 
# 	- In terms of proportion of artifactual SNVs within a sample, the max was 25%, mean was 0.37%
# 
# - Using a relaxed fp-cut of 5e-01: 
# 	- 5/204 samples (2.45%) was observed to possess artifacts. 
# 	- In terms of proportion of artifactual SNVs within a sample, the max was 25%, mean was 0.20%
# 
# **2025/10/27**
# 
# Out of 425 samples analyzed so far.
# 
# - Using a stringent fp-cut of 1e-08: 
# 	- 36/425 samples (8.5%) was observed to possess artifacts. 
# 	- In terms of proportion of artifactual SNVs within a sample, the max was 25%; mean was 0.34%
# 
# - Using a relaxed fp-cut of 5e-01: 
# 	- 13/425 samples (3.1%) was observed to possess artifacts. 
# 	- In terms of proportion of artifactual SNVs within a sample, the max was 25%; mean was 0.16%

# %% [markdown]
# ## FFPE-XML

# %%
get_res("xml", "ffpe", 5e-01)

# %%
get_res("xml", "ffpe", 5e-01, exclude_vus=True)

# %%
get_res("xml", "ffpe", 1e-08)

# %%
get_res("xml", "ffpe", 1e-08, exclude_vus=True)

# %% [markdown]
# ## OxoG-XML

# %%
get_res("xml", "oxog", 5e-01)

# %%
get_res("xml", "oxog", 5e-01, exclude_vus=True)

# %%
get_res("xml", "oxog", 1e-08)

# %%
get_res("xml", "oxog", 1e-08, exclude_vus=True)

# %% [markdown]
# #### Observations - XML variants
# 
# **2025/11/04**
# 
# 	FP-cut: 0.5
# 	Damage Type: FFPE
# 	Variant Source: XML
# 	Samples analyzed: 452
# 	Samples with ≥1 predicted artifact: 378 (83.6%)
# 	Mean proportion of artifactual SNVs: 16.87%
# 	Max proportion of artifactual SNVs: 62.50%
# 	Mean proportion of artifactual SNVs within samples with ≥1 predicted artifact: 20.17%
# 
# Results stratified by test type:
# 
# 	FoundationOneLiquidDx samples: 165
# 	FoundationOneLiquidDx max artifact proportion within samples: 42.86%
# 	FoundationOneLiquidDx mean artifact proportion within samples: 10.44%
# 	FoundationOneLiquidDx samples with >1 artifact: 119
# 	FoundationOneLiquidDx mean artifact proportion within samples with ≥1 detected artifacts: 14.47%
# 
# 	FoundationOne CDX samples: 287
# 	FoundationOne CDX max artifact proportion within samples: 62.50%
# 	FoundationOne CDX mean artifact proportion within samples: 20.57%
# 	FoundationOne CDX samples with >1 artifact: 259
# 	FoundationOne CDX mean artifact proportion within samples with ≥1 detected artifacts: 22.80%
# 
# 
# ------------------------
# 
# 	FP-cut: 1e-08
# 	Damage Type: FFPE
# 	Variant Source: XML
# 	Samples analyzed: 452
# 	Samples with ≥1 predicted artifact: 405 (89.6%)
# 	Mean proportion of artifactual SNVs: 20.34%
# 	Max proportion of artifactual SNVs: 71.43%
# 	Mean proportion of artifactual SNVs within samples with ≥1 predicted artifact: 22.70%
# 
# Results stratified by test type:
# 
# 	FoundationOneLiquidDx samples: 165
# 	FoundationOneLiquidDx max artifact proportion within samples: 50.00%
# 	FoundationOneLiquidDx mean artifact proportion within samples: 15.58%
# 	FoundationOneLiquidDx samples with >1 artifact: 139
# 	FoundationOneLiquidDx mean artifact proportion within samples with ≥1 detected artifacts: 18.50%
# 
# 	FoundationOne CDX samples: 287
# 	FoundationOne CDX max artifact proportion within samples: 71.43%
# 	FoundationOne CDX mean artifact proportion within samples: 23.07%
# 	FoundationOne CDX samples with >1 artifact: 266
# 	FoundationOne CDX mean artifact proportion within samples with ≥1 detected artifacts: 24.89%
# 
# 
# ----------------------
# 
# 	FP-cut: 0.5
# 	Damage Type: OXOG
# 	Variant Source: XML
# 	Samples analyzed: 452
# 	Samples with ≥1 predicted artifact: 191 (42.3%)
# 	Mean proportion of artifactual SNVs: 5.50%
# 	Max proportion of artifactual SNVs: 50.00%
# 	Mean proportion of artifactual SNVs within samples with ≥1 predicted artifact: 13.02%
# 
# Results stratified by test type:
# 
# 	FoundationOneLiquidDx samples: 165
# 	FoundationOneLiquidDx max artifact proportion within samples: 22.22%
# 	FoundationOneLiquidDx mean artifact proportion within samples: 3.18%
# 	FoundationOneLiquidDx samples with >1 artifact: 52
# 	FoundationOneLiquidDx mean artifact proportion within samples with ≥1 detected artifacts: 10.08%
# 
# 	FoundationOne CDX samples: 287
# 	FoundationOne CDX max artifact proportion within samples: 50.00%
# 	FoundationOne CDX mean artifact proportion within samples: 6.84%
# 	FoundationOne CDX samples with >1 artifact: 139
# 	FoundationOne CDX mean artifact proportion within samples with ≥1 detected artifacts: 14.12%
# 
# ---------------------
# 
# 	FP-cut: 1e-08
# 	Damage Type: OXOG
# 	Variant Source: XML
# 	Samples analyzed: 452
# 	Samples with ≥1 predicted artifact: 227 (50.2%)
# 	Mean proportion of artifactual SNVs: 7.08%
# 	Max proportion of artifactual SNVs: 50.00%
# 	Mean proportion of artifactual SNVs within samples with ≥1 predicted artifact: 14.10%
# 
# Results stratified by test type:
# 
# 	FoundationOneLiquidDx samples: 165
# 	FoundationOneLiquidDx max artifact proportion within samples: 25.00%
# 	FoundationOneLiquidDx mean artifact proportion within samples: 4.56%
# 	FoundationOneLiquidDx samples with >1 artifact: 68
# 	FoundationOneLiquidDx mean artifact proportion within samples with ≥1 detected artifacts: 11.06%
# 
# 	FoundationOne CDX samples: 287
# 	FoundationOne CDX max artifact proportion within samples: 50.00%
# 	FoundationOne CDX mean artifact proportion within samples: 8.53%
# 	FoundationOne CDX samples with >1 artifact: 159
# 	FoundationOne CDX mean artifact proportion within samples with ≥1 detected artifacts: 15.40%


