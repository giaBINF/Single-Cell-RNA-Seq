#!/bin/bash
#SBATCH --job-name=scrnaseq_partC
#SBATCH --output=partC_%j.out
#SBATCH --error=partC_%j.err
#SBATCH --time=00:30:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4

module load r

echo "========== PART C: Feature Plots =========="

Rscript --no-save --no-restore - <<'RSCRIPT'

library(Seurat)
library(ggplot2)

dir.create("plots", showWarnings = FALSE)

cat("Loading checkpoint from Part B...\n")
seurat_obj <- readRDS("checkpoints/after_partB.rds")
cat("Loaded:", ncol(seurat_obj), "cells\n")

pdf("plots/23_features_ISGs.pdf", width = 12, height = 10)
FeaturePlot(seurat_obj, features = c("Isg15", "Ifit1", "Oasl2", "Rsad2"))
dev.off()

pdf("plots/24_features_cytokines.pdf", width = 12, height = 10)
FeaturePlot(seurat_obj, features = c("Cxcl10", "Ccl2", "Il6", "Tnf"))
dev.off()

pdf("plots/25_features_antiviral.pdf", width = 12, height = 10)
FeaturePlot(seurat_obj, features = c("Ifnb1", "Ifng", "Stat1"))
dev.off()

cat("===== PART C COMPLETE =====\n")

RSCRIPT

echo "========== PART C DONE =========="
