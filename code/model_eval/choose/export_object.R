# CHOOSE (Zenodo 7083558, CHOOSE_ASD_Modules.rds = asd_chong_annot_withguides_v3.3modules_srt.rds of the authors'
# scripts): export raw RNA counts as dgCMatrix vectors (binary int32 / float32), gene and cell names, and the cell
# metadata, so that the clone test runs in Python.
# Usage: Rscript export_object.R <rds> <out dir>
args <- commandArgs(trailingOnly = TRUE)
suppressMessages({library(Seurat); library(Matrix)})
o <- readRDS(args[1]); d <- args[2]; dir.create(d, showWarnings = FALSE)
print(o); print(Assays(o))
m <- GetAssayData(o, assay = "RNA", slot = "counts")
m <- as(m, "dgCMatrix")
cat("counts", dim(m), length(m@x), "\n")
writeBin(as.integer(m@i), file.path(d, "i.bin"), size = 4)
writeBin(as.integer(m@p), file.path(d, "p.bin"), size = 4)
writeBin(as.numeric(m@x), file.path(d, "x.bin"), size = 4)
writeLines(rownames(m), file.path(d, "genes.txt")); writeLines(colnames(m), file.path(d, "cells.txt"))
md <- o@meta.data; md$cell <- rownames(md)
write.table(md, file.path(d, "meta.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
for (a in Assays(o)) cat("assay", a, dim(o[[a]]), "\n")
cat("meta columns:", colnames(md), "\n"); cat("DONE\n")
