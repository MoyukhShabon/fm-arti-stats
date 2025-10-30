import polars as pl
import os
import glob

def calc_arti_prop(artifact_pred_paths: list) -> pl.DataFrame:
	arti_prop = []

	for path in artifact_pred_paths:
		
		artifact_pred = pl.read_csv(path, separator="\t")

		snvs = artifact_pred.shape[0]
		artifacts = artifact_pred.filter(pl.col("pred") == False).shape[0]

		sample_stats = {
			"sample_name" : os.path.basename(path).split(".")[0],
			"n_snv" : snvs,
			"n_real" : (snvs - artifacts),
			"n_artifacts" : artifacts,
			"proportion" : (artifacts / snvs)
		}

		arti_prop.append(sample_stats)

	arti_prop = pl.DataFrame(arti_prop).sort("proportion", descending=True)
	return arti_prop

def get_stats(arti_prop: pl.DataFrame) -> None:
	mean_pct = arti_prop["proportion"].mean() * 100
	max_pct = arti_prop["proportion"].max() * 100
	n_samples = arti_prop.shape[0]
	n_with_artifacts = arti_prop.filter(pl.col("proportion") > 0).shape[0]

	print(f"Samples analyzed: {n_samples}")
	print(f"Samples with ≥1 predicted artifact: {n_with_artifacts} ({n_with_artifacts / n_samples * 100:.1f}%)")
	print(f"Mean proportion of artifactual SNVs: {mean_pct:.2f}%")
	print(f"Max proportion of artifactual SNVs: {max_pct:.2f}%")
	print("")

def get_res(damage_type: str, fp_cut: float) -> pl.DataFrame:
	
	search_patern = f"../{damage_type.lower()}-snvf/*/*.pred_fp-cut_{fp_cut:.0e}.tsv"
	arti_pred_paths = glob.glob(search_patern)
 
	if not arti_pred_paths:
		raise FileNotFoundError (f"No predictions made with fp-cut: {fp_cut:.0e} at {search_patern}")

	proportions = calc_arti_prop(arti_pred_paths)

	print(f"fp-cut: {fp_cut}")
	get_stats(proportions)

	proportions.write_csv(f"{damage_type}_proportions_per_sample.fp-cut_{fp_cut:.0e}.tsv", separator="\t")
	return proportions

# %% [markdown]
# ## FFPE

get_res("ffpe", 1e-08)

get_res("ffpe", 5e-01)

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
# ## OxoG

get_res("oxog", 1e-08)

get_res("oxog", 5e-01)

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


