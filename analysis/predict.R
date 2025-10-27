#!/usr/bin/env Rscript

library(io)
library(argparser)
source("../../common-ffpe-snvf/R/eval.R")

p <- arg_parser("Predict Artifactual SNVs via a FDR adjusted False Positive threshold")
p <- add_argument(p, "--fp.cut", help = "False positive cutoff (numeric, between 0 and 1)", default = 0.5)
argv <- parse_args(p)
fp.cut <- as.numeric(argv$fp.cut)
if (is.na(fp.cut) || fp.cut < 0 || fp.cut > 1) {
    stop("Invalid --fp.cut; must be a number between 0 and 1")
}

dir.create("../annot", showWarnings=FALSE)

snvf_paths <- Sys.glob("../*/*/*.mobsnvf.*.snv")

# Get file info for all paths at once
file_sizes <- file.info(snvf_paths)$size

# Separate valid and empty files
empty_files <- snvf_paths[is.na(file_sizes) | file_sizes == 0]
valid_files <- snvf_paths[!is.na(file_sizes) & file_sizes > 0]

# Save empty files list if any
if (length(empty_files) > 0) {
    writeLines(empty_files, "../annot/no_variants.txt")
    message(sprintf("Found %d empty/missing files", length(empty_files)))
}

for (path in valid_files){

    message(sprintf("Processing: %s", basename(path)))

    out_path <- gsub("\\.snv", sprintf("\\.pred_fp-cut_%s.tsv", format(fp.cut, scientific = TRUE)), path)
    snv_pred <- fdr_cut_pred(read.delim(path), "FOBP", fp.cut)

    qwrite(snv_pred, out_path)
    message(sprintf("Prediction for %s written to: %s", basename(path), out_path))

}
