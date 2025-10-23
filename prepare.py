#!/usr/bin/env python
import os
import glob
import polars as pl

# %%
def assign_if_exist(path: str):
	if not os.path.exists(path):
		raise FileNotFoundError(f"File not found: {path}")
	else:
		return path

# %%
ffpe_outdir = "ffpe-snvf"
oxog_outdir = "oxog-snvf"

# %%
ref_path = os.path.abspath(assign_if_exist("ref/hg19/ucsc.hg19.fasta"))

# %%
bam_paths = glob.glob("data/*.bam")
vcf_paths = glob.glob("data/*.vcf")


bam_table = pl.DataFrame({
	"sample_name" : [os.path.basename(path).split("*")[0].split("_")[0] for path in bam_paths],	
	"bam_path": [os.path.abspath(path) for path in bam_paths]
})


vcf_table = pl.DataFrame({
	"sample_name" : [os.path.basename(path).removesuffix(".vcf") for path in vcf_paths],	
	"vcf_path": [os.path.abspath(path) for path in vcf_paths]
})


bam_vcf_table = bam_table.join(vcf_table, on="sample_name", how = "inner")
bam_vcf_table

# %%
templates = ["ffpe-snvf/mobsnvf.ffpe.sh.template", "oxog-snvf/mobsnvf.oxog.sh.template"]

for i, sample_name in enumerate(bam_vcf_table["sample_name"]):
	print(f"Creating scripts for {sample_name}")

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
