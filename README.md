# FM-ARTI-STATS

Foundation Medicine Artifact Status

## Objective

Here we aim to analyze the extent of artifacts - FFPE, OxoG etc that may be present in foundation medicine dataset.

## Method

MOBSNVF is used to determine the extent of FFPE and OxoG artifacts under the assumption that it has sufficient precision i.e the artifacts identified are True Positive.

Analysis is performed using variants from both VCFs and XMLs for each sample with their corresponding BAM.

### XML 

SNVs are parsed from the XML into an SNV format before being passed on to each sample

### VCF

Raw VCFs are used as input into MOBSNVF

### Prediction

Raw results from MOBSNVF are stratified to predict true variants. 

The scores from MOBSNVF is first FDR adjusted and then these scores are used for making predictions based on a False Positive Cutoff threshold.

Two thresholds of False Positive Cutoff (FP-cut) is used in the analysis as of now:
- 5e-01
- 1e-08

### Artifact Proportion

These predictions are then used to calculate the proportion of artifacts in each sample and the overall statistics.


## Observations

The BAM headers states that read are aligned against the hg19 reference. Therefore the UCSC hg19 reference was downloaded from the broad institute.

```
gs://gatk-legacy-bundles
```


## Replication

**Note:** All directory paths mentioned below are relative to the repository root.

1. The repository works under the assumption that the Foundation Medicine Data i.e BAM and VCFs placed in the `data` directory.

    A `link.sh` script is placed in this directory which can be modified to create a symbolic link of the data from the path where it resides.

2. Download the reference genome from the Broad Institute's legacy bundles at `gs://gatk-legacy-bundles`. 

Navigate to the `ref` directory and run:

```bash
bash get.sh
```

3. Navigate to the repository root and run `prepare-vcf.py` and `prepare-xml.py`


```bash
python prepare-vcf.py
python prepare-xml.py
```

This will create batch scripts for FFPE and OxoG artifact filtering in the `vcf-ffpe-snvf` and `vcf-oxog-snvf` directories respectively based on the BAM and VCF data available in the `data` directory.

Similarly, batch scripts will also be created in the `xml-ffpe-snvf` and `xml-oxog-snvf` directories based on the BAMs and SNV information contained within the XMLs.

The variant data parsed from the XML are saved in a tabular format in the `xml-snvs`. The a `<sample_name>.tsv` and `<sample_name.snv>` file is created based on for each sample's representative XML. The SNV file only contains variant info i.e CHROM, POS, REF, and ALT columns. The TSV file includes additional information such as Allele Fraction, Gene, Protein Effect, Strand, Functional Effect etc.

4. Navigate to `vcf-ffpe-snvf` directory and run `filter.sh` to perform MOSBNVF filtering for FFPE artifacts on the variants present in the VCFs.

```bash
    bash filter.sh
```

5. Navigate to the `vcf-oxog-snvf` directory and run `filter.sh` to perform MOSBNVF filtering for OxoG artifacts on the variants present in the VCFs.

```bash
    bash filter.sh
```

6. Navigate to `xml-ffpe-snvf` directory and run `filter.sh` to perform MOSBNVF filtering for FFPE artifacts on the variants present in the XMLs.

```bash
    bash filter.sh
```

7. Navigate to the `xml-oxog-snvf` directory and run `filter.sh` to perform MOSBNVF filtering for OxoG artifacts on the variants present in the XMLs.

```bash
    bash filter.sh
```

8. Navigate to the `analysis` directory and run `predict.R` with specific values of fp-cut.

```
Rscript predict.R --fp.cut 0.5
Rscript predict.R --fp.cut 1e-08
```

This creates predictions based on the specified False Positive Cut threshold applied after FDR correcting the scores of MOBSNVF. The predictions are saved to the same directories as the MOBSNVF results mentioned in step 4-7.


9. Run the `proportions.py` under the same analysis directory:

This calculates the proportion of artifacts in each sample and gives an overall summary. The table with the proportions are saved in the same directory. E.g. `ffpe_proportions_per_sample.xml.fp-cut_1e-08.tsv`, `oxog_proportions_per_sample.vcf.fp-cut_5e-01.tsv` etc.

## Issues

It was seen that some VCFs have multiple BAMs. It was seen that the <sample_name>_DNA.bam is significantly smaller than the <sample_name>*US\<nnnnnnn>.sorted.bam. Sumedha was contacted for elaboration of the data. As of right now the larger bam is being used for analysis.

![multiple_bams_per_sample](.md_images/multiple_bam_per_sample.png)

Some VCFs do not have associated BAMs

![samples_no_bams](.md_images/no_bam_samples.png)

