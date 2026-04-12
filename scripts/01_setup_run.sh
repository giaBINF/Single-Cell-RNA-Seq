#!/bin/bash
#SBATCH --job-name=scrnaseq_partA
#SBATCH --output=partA_%j.out
#SBATCH --error=partA_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=4

module load r

echo "========== PART A: Load, QC, Normalize, Cluster =========="

Rscript --no-save --no-restore - <<'RSCRIPT'

library(Seurat)
library(ggplot2)

dir.create("plots", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("checkpoints", showWarnings = FALSE)

cat("Loading Seurat object...\n")
seurat_obj <- readRDS("seurat_ass4.rds")
cat("Loaded:", ncol(seurat_obj), "cells x", nrow(seurat_obj), "genes\n")

cat("\n--- Metadata columns ---\n")
print(colnames(seurat_obj@meta.data))
cat("\n--- Time points ---\n")
print(table(seurat_obj$time))
cat("\n--- Tissue types ---\n")
print(table(seurat_obj$organ_custom))
cat("\n--- Design matrix ---\n")
print(table(seurat_obj$organ_custom, seurat_obj$time))
cat("\n--- Missing mouse_id count ---\n")
print(sum(is.na(seurat_obj$mouse_id) | seurat_obj$mouse_id == ""))

cat("\nRunning QC...\n")
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")

pdf("plots/01_qc_violin.pdf", width = 12, height = 5)
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3, pt.size = 0)
dev.off()

pdf("plots/02_qc_scatter_mt.pdf", width = 8, height = 6)
FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
dev.off()

pdf("plots/03_qc_scatter_features.pdf", width = 8, height = 6)
FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

seurat_obj <- subset(seurat_obj,
                     subset = nFeature_RNA > 200 &
                              nFeature_RNA < 6000 &
                              percent.mt < 20)
cat("After QC:", ncol(seurat_obj), "cells remain\n")

cat("\nNormalizing...\n")
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize")
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst",
                                   nfeatures = 2000)
cat("Scaling variable features...\n")
seurat_obj <- ScaleData(seurat_obj)

cat("Running PCA...\n")
seurat_obj <- RunPCA(seurat_obj,
                     features = VariableFeatures(object = seurat_obj))

pdf("plots/04_elbow_plot.pdf", width = 8, height = 6)
ElbowPlot(seurat_obj, ndims = 40)
dev.off()

n_pcs <- 25
cat("Clustering (resolution 0.3)...\n")
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:n_pcs)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.3)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:n_pcs)
cat("Number of clusters:", length(levels(Idents(seurat_obj))), "\n")

pdf("plots/05_umap_clusters.pdf", width = 10, height = 8)
DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = 0.3) +
  ggtitle("UMAP — All Clusters")
dev.off()

pdf("plots/06_umap_by_tissue.pdf", width = 18, height = 6)
DimPlot(seurat_obj, reduction = "umap", split.by = "organ_custom", label = TRUE)
dev.off()

pdf("plots/07_umap_by_time.pdf", width = 24, height = 6)
DimPlot(seurat_obj, reduction = "umap", split.by = "time", label = TRUE)
dev.off()

pdf("plots/08_umap_batch_origident.pdf", width = 10, height = 8)
DimPlot(seurat_obj, reduction = "umap", group.by = "orig.ident")
dev.off()

pdf("plots/09_umap_batch_organ.pdf", width = 10, height = 8)
DimPlot(seurat_obj, reduction = "umap", group.by = "organ_custom")
dev.off()

cat("Saving checkpoint...\n")
saveRDS(seurat_obj, "checkpoints/after_partA.rds")
cat("===== PART A COMPLETE =====\n")

RSCRIPT

echo "========== PART A DONE =========="
