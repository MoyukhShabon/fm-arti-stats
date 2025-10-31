#!/usr/bin/env python
import polars as pl
import os
import glob
from lxml import etree

def get_test_type(xml_path: int) -> str:
	
	if not os.path.exists(xml_path):
		raise FileNotFoundError(f"No file found at {xml_path}")
	
	tree = etree.parse(xml_path)
	root = tree.getroot()
	return root.xpath("//TestType")[0].text


def calc_arti_prop(artifact_pred_paths: list) -> pl.DataFrame:
	arti_prop = []

	for path in artifact_pred_paths:
		
		artifact_pred = pl.read_csv(path, separator="\t")

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

def get_stats(arti_prop: pl.DataFrame, damage_type: str = None) -> None:
	
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

	print(f"Samples analyzed: {n_samples}")
	print(f"Samples with ≥1 predicted artifact: {n_with_artifacts} ({n_with_artifacts / n_samples * 100:.1f}%)")
	print(f"Mean proportion of artifactual SNVs: {mean_pct:.2f}%")
	print(f"Max proportion of artifactual SNVs: {max_pct:.2f}%")
	print(f"Mean proportion of artifactual SNVs within samples with ≥1 predicted artifact: {mean_pct_arti_samples:.2f}%")
 
	print("Results stratified by test type:")
	for i in range(prop_test_type.shape[0]):

		test_type = prop_test_type[i, "test_type"]

		print(f"\t{test_type} samples: {prop_test_type[i, "count"]}")
		print(f"\t{test_type} max artifact proportion within samples: {prop_test_type[i, "max_proportion"] * 100:.2f}%")
		print(f"\t{test_type} mean artifact proportion within samples: {prop_test_type[i, "mean_proportion"] * 100:.2f}%")

		prop_test_type_with_artifact_filtered = prop_test_type_with_artifact.filter(pl.col("test_type") == test_type)

		if prop_test_type_with_artifact_filtered.shape[0]:
			print(f"\t{test_type} samples with >1 artifact: {prop_test_type_with_artifact_filtered[0, "count_>=1_artifact"]}")
			print(f"\t{test_type} mean artifact proportion within samples with ≥1 detected artifacts: {prop_test_type_with_artifact_filtered[0, "mean_proportion"] * 100:.2f}%")

		## Create spacing
		print()
	print()

def get_res(damage_type: str, fp_cut: float) -> pl.DataFrame:
	
	search_pattern = f"../{damage_type.lower()}-snvf/*/*.pred_fp-cut_{fp_cut:.0e}.tsv"
	arti_pred_paths = glob.glob(search_pattern)
 
	if not arti_pred_paths:
		raise FileNotFoundError (f"No predictions made with fp-cut: {fp_cut:.0e} at {search_pattern}")

	proportions = calc_arti_prop(arti_pred_paths)

	print(f"FP-cut: {fp_cut}")
	get_stats(proportions, damage_type=damage_type)

	proportions.write_csv(f"{damage_type}_proportions_per_sample.fp-cut_{fp_cut:.0e}.tsv", separator="\t")
	return proportions

# ## FFPE
get_res("ffpe", 1e-08)

get_res("ffpe", 5e-01)


# ## OxoG
get_res("oxog", 1e-08)

get_res("oxog", 5e-01)
