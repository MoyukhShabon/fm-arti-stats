#!/usr/bin/env bash

import polars as pl
from lxml import etree
import os
import glob
import re

# %% [markdown]
# ## Functions

# %%
def parse_snvs_from_xml(xml_path: str) -> pl.DataFrame:
    """
    Parses an XML report to extract Single Nucleotide Variants (SNVs) and
    returns them as a Polars DataFrame.

    This function identifies SNVs by parsing the 'cds-effect' attribute 
    searching for a single nucleotide change (e.g., '829G>A').
    This excludes INDELs (insertions/deletions) and other complex variants.

    Args:
        xml_path: The file path to the input XML report.

    Returns:
        A Polars DataFrame containing the parsed SNVs with columns:
        'chrom', 'pos', 'ref', 'alt', 'gene', 'depth', 'protein_effect',
        'allele_fraction', 'functional_effect', and 'transcript'.
        Returns an empty DataFrame if no SNVs are found.
    """

	# ---- Setup ---

    if not os.path.exists(xml_path):
        raise FileNotFoundError(f"No file found at {xml_path}")

    try:
        tree = etree.parse(xml_path)
        root = tree.getroot()
    except etree.XMLSyntaxError as e:
        raise ValueError(f"Error parsing XML file: {e}")


    # The <variant-report> and its children use a default namespace.
    # We must define a prefix (e.g., 'vr') to use in our XPath queries.
    namespaces = {
        'vr': 'http://foundationmedicine.com/compbio/variant-report-external'
    }

	# ---- SNV extraction -----
    variants_data = []

    # Regex is used to identify SNVs from a cds-effect string within the XML (e.g., cds_effect"829G>A").
    # It captures a single reference nucleotide [ACGTN] followed by '>' and a single alternate nucleotide [ACGTN]. 
    # This inherently filters out 'ins' and 'del'.
    # Note: lxml automatically handles XML entities, so '829G&gt;A' is read as '829G>A'.
    snv_regex = re.compile(r'^[\d\+\-\*\_]+([ACGTN])>([ACGTN])$')

    # Use XPath to find all 'short-variant' elements within the defined namespace.
    for variant in root.xpath('.//vr:short-variant', namespaces=namespaces):
        cds_effect_str = variant.get('cds-effect')
        position_str = variant.get('position')
        equivocal_str = variant.get('equivocal')

        # Ensure essential attributes are present before proceeding
        if not all([position_str, cds_effect_str]):
            continue

        # --- SNV Filtering and Parsing ---
        # Use the regex to check if the cds-effect describes a simple SNV.
        match = snv_regex.match(cds_effect_str)
        
        # If it's a match, we've found an SNV.
        if match:
            try:
                # The regex match gives us REF and ALT directly.
                ref, alt = match.groups()

                # Parse chromosome and position
                chrom, pos_str = position_str.split(':')
                pos = int(pos_str)
                is_equivocal = equivocal_str.lower() == 'true'

                # Extract other relevant information, converting types where necessary
                variant_dict = {
                    'chrom': chrom,
                    'pos': pos,
                    'ref': ref,
                    'alt': alt,
                    'gene': variant.get('gene'),
                    'depth': int(variant.get('depth')),
                    'protein_effect': variant.get('protein-effect'),
                    'allele_fraction': float(variant.get('allele-fraction')),
                    'functional_effect': variant.get('functional-effect'),
                    'transcript': variant.get('transcript'),
                    'strand': variant.get('strand'),
                    'equivocal': is_equivocal,
                    # 'test_type': root.xpath("//TestType")[0].text
                }
                variants_data.append(variant_dict)

            except (ValueError, IndexError):
                # Skip this variant if parsing fails (e.g., malformed position)
                print(f"Skipping variant due to parsing error: {etree.tostring(variant)}")
                continue
    

    if not variants_data:
        print(f"No SNVs detected: {xml_path}")
        return pl.DataFrame() # Return an empty DataFrame if no SNVs were found


    df = pl.DataFrame(variants_data)
    
    # ------ formatting -------
    
    chrom_sort_key = (
		pl.when(pl.col("chrom") == "X").then(pl.lit(23, dtype=pl.Int64))
		.when(pl.col("chrom") == "Y").then(pl.lit(24, dtype=pl.Int64))
		.when((pl.col("chrom") == "M") | (pl.col("chrom") == "MT")).then(pl.lit(25, dtype=pl.Int64))
		.otherwise(pl.col("chrom").str.extract(r"(\d+)", 1).cast(pl.Int64))
		.fill_null(99) # Place any other contigs at the very end
	)
    
    # Ensure column order is consistent
    schema = [
        'chrom', 'pos', 'ref', 'alt', 
        'gene', 'depth', 'protein_effect',
        'allele_fraction', 'functional_effect', 
        'transcript', 'strand', 'equivocal', 
        # "test_type",
    ]
    df = df.select(schema).sort(chrom_sort_key, pl.col("pos"))
    
    return df


# %%
def assign_if_exist(path: str) -> str:
	if not os.path.exists(path):
		raise FileNotFoundError(f"File not found: {path}")
	else:
		return path

# %% [markdown]
# ## Parse XML to SNV

# %%
outdir_root = "xml-snvs"

xml_paths = sorted(glob.glob("data/*.xml"))

# %%
no_variants = []

for i, path in enumerate(xml_paths):
	
	print(f"{i+1}. Processing {path}")	
 
	sample_name = os.path.basename(path).replace(".xml", "")
	outdir = f"{outdir_root}/{sample_name}"
	os.makedirs(outdir, exist_ok=True)
 
	snvs = parse_snvs_from_xml(path)
	
	if snvs.is_empty():
		print(f'Dataframe for "{sample_name}" is empty. Skipping writing to disk')
		no_variants.append(path)
		continue
 
	snvs.write_csv(f"{outdir}/{sample_name}.tsv", separator="\t")
	snvs.select(["chrom", "pos", "ref", "alt"]).write_csv(f"{outdir}/{sample_name}.snv", separator="\t")
	

pl.DataFrame(no_variants).write_csv("annot/xml-no_snv.txt", include_header=False, separator="\t")

# %% [markdown]
# ## Make FFPE-SNVF scripts

# %%
ffpe_outdir = "xml-ffpe-snvf"
oxog_outdir = "xml-oxog-snvf"

ref_path = os.path.abspath(assign_if_exist("ref/hg19/ucsc.hg19.fasta"))

bam_paths = glob.glob("data/*.bam")
snv_paths = glob.glob("xml-snvs/*/*.snv")


bam_table = pl.DataFrame({
	"sample_name" : [os.path.basename(path).split("*")[0].split("_")[0] for path in bam_paths],	
	"bam_path": [os.path.abspath(path) for path in bam_paths]
})


snv_table = pl.DataFrame({
	"sample_name" : [os.path.basename(path).removesuffix(".snv") for path in snv_paths],	
	"snv_path": [os.path.abspath(path) for path in snv_paths]
})

no_bam = snv_table.join(bam_table, on="sample_name", how = "anti")
no_bam.write_csv("annot/xml-no_bam.tsv", separator="\t")

bam_snv_table = bam_table.join(snv_table, on="sample_name", how = "inner")
# In cases where a sample has two BAMs, only keep the BAM which is <sample_id>*US<nnnnnnn>.sorted.bam
bam_snv_table = bam_snv_table.filter(~(bam_snv_table["sample_name"].is_duplicated() & pl.col("bam_path").str.contains("_DNA.bam")))
bam_snv_table


# %%

templates = ["xml-ffpe-snvf/mobsnvf.ffpe.sh.template", "xml-oxog-snvf/mobsnvf.oxog.sh.template"]

for i, sample_name in enumerate(bam_snv_table["sample_name"]):
	print(f"Creating scripts for {sample_name}")

	bam_path = bam_snv_table[i, "bam_path"]
	snv_path = bam_snv_table[i, "snv_path"]

	for template in templates:
	
		filtered_outdir = os.path.abspath(f"{template.split("/")[0]}")
		script_outdir = f"{filtered_outdir}/scripts"
		os.makedirs(script_outdir, exist_ok=True)
		
		with open(template, "r") as t:
			script = t.read()

		new_script = script.replace("$1", bam_path).replace("$2", snv_path).replace("$3", ref_path).replace("$4", f"{filtered_outdir}/{sample_name}")

		script_outpath = f"{script_outdir}/{sample_name}.sh"

		with open(script_outpath, "w") as f:
			f.write(new_script)


