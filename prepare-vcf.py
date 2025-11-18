#!/usr/bin/env python
import os
import glob
import polars as pl

# %%
def return_path_if_exists(path: str) -> str:
	if not os.path.exists(path):
		raise FileNotFoundError(f"File not found: {path}")
	else:
		return path

# %%
ffpe_outdir = "vcf-ffpe-snvf"
oxog_outdir = "vcf-oxog-snvf"

# %%
ref_path = os.path.abspath(return_path_if_exists("ref/hg19/ucsc.hg19.fasta"))

# %%
## Get BAM and VCF paths
bam_paths = sorted(glob.glob("data/*.bam"))
vcf_paths = sorted(glob.glob("data/*.vcf"))
print(f"Found {len(bam_paths)} BAM files and {len(vcf_paths)} VCF files")

## Make BAM and VCF tables with path and sample name inferred from file names
bam_table = pl.DataFrame({
	"sample_name" : [os.path.basename(path).split("*")[0].split("_")[0] for path in bam_paths],	
	"bam_path": [os.path.abspath(path) for path in bam_paths]
})


vcf_table = pl.DataFrame({
	"sample_name" : [os.path.basename(path).removesuffix(".vcf") for path in vcf_paths],	
	"vcf_path": [os.path.abspath(path) for path in vcf_paths]
})

## Go through each VCF and check if there are any variants and create a table of VCFs with no variants
print("Checking if there are VCFs without variants...")
vcf_no_variants = []
for i, sample_name in enumerate(vcf_table["sample_name"]):
    print(f"{i+1}. {sample_name}")
    df = pl.read_csv(vcf_table[i, "vcf_path"], separator="\t", comment_prefix="##")
    
    if df.is_empty():
        print(f"\t{sample_name} has no variants")
        vcf_no_variants.append(vcf_table[i])
    
vcf_no_variants = pl.concat(vcf_no_variants)
vcf_no_variants.write_csv("annot/vcf-no_variants.tsv", separator = "\t")
print(f"Found {vcf_no_variants.height} VCFs with no variants \nTable written to annot/vcf-no_variants.tsv \n")


## Find VCFs with no corresponding BAMs
no_bam = vcf_table.join(bam_table, on="sample_name", how = "anti")
no_bam.write_csv("annot/vcf-no_bam.tsv", separator = "\t")

## Join BAM and VCF tables to get only samples with both BAM and VCF
bam_vcf_table = bam_table.join(vcf_table, on="sample_name", how = "inner")

# In cases where a sample has two BAMs, only keep the BAM which is <sample_id>*US<nnnnnnn>.sorted.bam
bam_vcf_table = bam_vcf_table.filter(~(bam_vcf_table["sample_name"].is_duplicated() & pl.col("bam_path").str.contains("_DNA.bam")))
## Remove samples with no variants from bam_vcf_table
bam_vcf_table = bam_vcf_table.join(vcf_no_variants, on="sample_name", how="anti")
bam_vcf_table

# %%
os.makedirs("annot", exist_ok=True)
bam_vcf_table.write_csv("annot/bam_vcf_path.absolute.tsv", separator="\t")

# %%
## list containing path to MOBSNVF OXOG and FFPE templates
templates = ["vcf-ffpe-snvf/mobsnvf.ffpe.sh.template", "vcf-oxog-snvf/mobsnvf.oxog.sh.template"]

## Create scripts for each sample based on the templates
for i, sample_name in enumerate(bam_vcf_table["sample_name"]):
	print(f"Creating scripts for {i+1}. {sample_name}")

	bam_path = bam_vcf_table[i, "bam_path"]
	vcf_path = bam_vcf_table[i, "vcf_path"]

	for template in templates:
	
		filtered_outdir = os.path.abspath(f"{template.split("/")[0]}")
		script_outdir = f"{filtered_outdir}/scripts"
		os.makedirs(script_outdir, exist_ok=True)
		
		with open(template, "r") as t:
			script = t.read()

		new_script = script.replace("$1", bam_path).replace("$2", vcf_path).replace("$3", ref_path).replace("$4", f"{filtered_outdir}/{sample_name}")

		script_outpath = f"{script_outdir}/{sample_name}.sh"

		with open(script_outpath, "w") as f:
			f.write(new_script)



