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

snvf_paths <- Sys.glob("../*/*/*.mobsnvf.*.snv")

for (path in snvf_paths){

    fi <- file.info(path)
    if (is.na(fi$size) || fi$size == 0) {
        message(sprintf("Skipping missing/empty file: %s", basename(path)))
        next
    }

    message(sprintf("Processing: %s", basename(path)))

    out_path <- gsub("\\.snv", sprintf("\\.pred_fp-cut_%s.tsv", format(fp.cut, scientific = TRUE)), path)
    snv_pred <- fdr_cut_pred(read.delim(path), "FOBP", fp.cut)

    qwrite(snv_pred, out_path)
    message(sprintf("Prediction for %s written to: %s", basename(path), out_path))

}


