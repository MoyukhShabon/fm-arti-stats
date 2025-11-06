#!/usr/bin/env python
import polars as pl

# %% [markdown]
# ### Functions

def add_suffix(name: str, suffix: str) -> str:
	exceptions = ["sample_name", "test_type"]
	if name in exceptions:
		return name
	else:
		return f"{name}{suffix}"

def load_proportion_summary(damage_type: str, fp_cut: float, source: str, vus: str) -> pl.DataFrame:
    return pl.read_csv(f"{damage_type}_proportions_per_sample.{source}{vus}.fp-cut_{fp_cut:.0e}.tsv", separator="\t").rename(lambda x : add_suffix(x, f"_{source}{vus}"))

# %% [markdown]
# ### Combine Tables

## FFPE FP Cut - 0.5
vcf = load_proportion_summary("ffpe", 0.5, "vcf", "")
xml = load_proportion_summary("ffpe", 0.5, "xml", "")
xml_no_vus = load_proportion_summary("ffpe", 0.5, "xml", ".no_vus")

compiled = (
    vcf
    .join(xml, on=["sample_name", "test_type"], how="left")
    .join(xml_no_vus, on=["sample_name", "test_type"], how="left")
)

compiled.write_csv("ffpe_combined.vcf.xml.xml-no_vus.5e-01.tsv", separator="\t")

## FFPE FP Cut - 1e-08

vcf = load_proportion_summary("ffpe", 1e-08, "vcf", "")
xml = load_proportion_summary("ffpe", 1e-08, "xml", "")
xml_no_vus = load_proportion_summary("ffpe", 1e-08, "xml", ".no_vus")

compiled = (
    vcf
    .join(xml, on=["sample_name", "test_type"], how="left")
    .join(xml_no_vus, on=["sample_name", "test_type"], how="left")
)

compiled.write_csv("ffpe_combined.vcf.xml.xml-no_vus.1e-08.tsv", separator="\t")

## OXOG FP Cut - 0.5


vcf = load_proportion_summary("oxog", 0.5, "vcf", "")
xml = load_proportion_summary("oxog", 0.5, "xml", "")
xml_no_vus = load_proportion_summary("oxog", 0.5, "xml", ".no_vus")

compiled = (
    vcf
    .join(xml, on=["sample_name", "test_type"], how="left")
    .join(xml_no_vus, on=["sample_name", "test_type"], how="left")
)

compiled.write_csv("oxog_combined.vcf.xml.xml-no_vus.5e-01.tsv", separator="\t")

## OXOG FP Cut - 1e-08

vcf = load_proportion_summary("oxog", 1e-08, "vcf", "")
xml = load_proportion_summary("oxog", 1e-08, "xml", "")
xml_no_vus = load_proportion_summary("oxog", 1e-08, "xml", ".no_vus")

compiled = (
    vcf
    .join(xml, on=["sample_name", "test_type"], how="left")
    .join(xml_no_vus, on=["sample_name", "test_type"], how="left")
)

compiled.write_csv("oxog_combined.vcf.xml.xml-no_vus.1e-08.tsv", separator="\t")


