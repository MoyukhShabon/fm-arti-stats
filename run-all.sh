#!/bin/bash

# set -euo pipefail

echo -e "linking Foundation medicine data to data/ directory \n"
cd data
bash link.sh
cd -

echo -e "\nPreparing VCF SNVF scripts \n"
python prepare-vcf.py

echo -e "\nRunning VCF FFPE SNVF\n"
cd vcf-ffpe-snvf/scripts
dlazy .
cd -

echo -e "\nRunning VCF OXOG SNVF\n"
cd vcf-oxog-snvf/scripts
dlazy .
cd -

echo -e "\nPreparing XML SNV Tables and making XML SNVF results by subsetting of VCF SNVF results \n"
python prepare-xml.py

echo -e "\nCreating VUS excluded set from FFPE SNVF results  \n"
python filter_vus.py


echo -e "\nPredicting Artifacts \n"
Rscript predict.R 
Rscript predict.R --fp.cut 1e-08


echo -e "\nCalculating Proportion of artifacts per sample \n"
cd analysis
python proportion.py | tee summary_statistics.txt

echo -e "\nCombining artifact proportion tables \n"
python compilation.py

echo -e "\nCollecting all artifacts into a table \n"
python collect-artifacts.py
python collect-artifacts.py --fp-cut 0.5

