#!/usr/bin/env Rscript

library(io)
source("../../common-ffpe-snvf/R/eval.R")

snvf_paths <- Sys.glob("../*/*/*.mobsnvf.*.snv")


for (path in snvf_paths){

    fi <- file.info(path)
    if (is.na(fi$size) || fi$size == 0) {
        message(sprintf("Skipping missing/empty file: %s", basename(path)))
        next
    }

    message(sprintf("Processing: %s", basename(path)))

    out_path <- gsub("\\.snv", "\\.pred_fp-cut_5e-1.tsv", path)
    snv_pred <- fdr_cut_pred(read.delim(path), "FOBP", fp.cut = 0.5)

    qwrite(snv_pred, out_path)
    message(sprintf("Prediction for %s written to: %s", basename(path), out_path))

}


