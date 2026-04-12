#!/bin/bash
#SBATCH --job-name=scrnaseq_partB
#SBATCH --output=partB_%j.out
#SBATCH --error=partB_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=4

module load r

echo "========== PART B: FindAllMarkers + Annotation =========="

Rscript --no-save --no-restore - <<'RSCRIPT'

library(Seurat)
library(ggplot2)
library(dplyr)
library(SingleR)
library(celldex)

dir.create("plots", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("checkpoints", showWarnings = FALSE)

cat("Loading checkpoint from Part A...\n")
seurat_obj <- readRDS("checkpoints/after_partA.rds")
cat("Loaded:", ncol(seurat_obj), "cells,",
    length(levels(Idents(seurat_obj))), "clusters\n")

# --- FindAllMarkers ---
cat("\nFinding all markers...\n")
all.markers <- FindAllMarkers(seurat_obj,
                              only.pos        = TRUE,
                              min.pct         = 0.25,
                              logfc.threshold = 0.25)

top5 <- all.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

write.csv(all.markers, "results/all_markers.csv", row.names = FALSE)
write.csv(top5, "results/top5_markers_per_cluster.csv", row.names = FALSE)

cat("Top 5 markers per cluster:\n")
print(top5, n = Inf)

pdf("plots/10_heatmap_top_markers.pdf", width = 16, height = 12)
DoHeatmap(seurat_obj, features = top5$gene, size = 3) + NoLegend()
dev.off()

# --- SingleR ---
cat("\nRunning SingleR annotation...\n")
ref <- celldex::ImmGenData()
expr_matrix <- GetAssayData(seurat_obj, layer = "data")

singler_results <- SingleR(test     = expr_matrix,
                           ref      = ref,
                           labels   = ref$label.main,
                           clusters = Idents(seurat_obj))

cat("SingleR labels:\n")
print(singler_results$labels)

label_map <- singler_results$labels
names(label_map) <- rownames(singler_results)
singler_label <- label_map[as.character(Idents(seurat_obj))]
seurat_obj$singler_label <- unname(singler_label)

pdf("plots/11_umap_singler.pdf", width = 12, height = 8)
DimPlot(seurat_obj, group.by = "singler_label", label = TRUE, repel = TRUE) +
  ggtitle("Automated Annotation (SingleR — ImmGen)")
dev.off()

# --- FeaturePlots ---
cat("\nGenerating annotation FeaturePlots...\n")

pdf("plots/12_features_tcells.pdf", width = 12, height = 10)
FeaturePlot(seurat_obj, features = c("Cd3d", "Cd3e", "Cd4", "Cd8a"))
dev.off()

pdf("plots/13_features_bcells.pdf", width = 12, height = 5)
FeaturePlot(seurat_obj, features = c("Cd79a", "Ms4a1"))
dev.off()

pdf("plots/14_features_macrophages.pdf", width = 12, height = 10)
FeaturePlot(seurat_obj, features = c("Adgre1", "Cd68", "Lyz2"))
dev.off()

pdf("plots/15_features_neutrophils.pdf", width = 12, height = 5)
FeaturePlot(seurat_obj, features = c("S100a8", "S100a9"))
dev.off()

pdf("plots/16_features_nk.pdf", width = 12, height = 5)
FeaturePlot(seurat_obj, features = c("Nkg7", "Gzma"))
dev.off()

pdf("plots/17_features_dc.pdf", width = 12, height = 10)
FeaturePlot(seurat_obj, features = c("Itgax", "Flt3", "H2-Ab1"))
dev.off()

pdf("plots/18_features_epithelial.pdf", width = 12, height = 10)
FeaturePlot(seurat_obj, features = c("Foxj1", "Scgb1a1", "Krt5", "Krt14"))
dev.off()

pdf("plots/19_features_olfactory.pdf", width = 12, height = 5)
FeaturePlot(seurat_obj, features = c("Omp", "Cyp2g1"))
dev.off()

pdf("plots/20_features_stromal.pdf", width = 12, height = 10)
FeaturePlot(seurat_obj, features = c("Col1a1", "Dcn", "Pecam1"))
dev.off()

pdf("plots/21_dotplot_markers.pdf", width = 16, height = 10)
DotPlot(seurat_obj,
        features = c("Cd3d", "Cd4", "Cd8a", "Cd79a", "Adgre1", "Lyz2",
                      "S100a8", "Nkg7", "Itgax", "H2-Ab1",
                      "Foxj1", "Scgb1a1", "Krt5", "Omp",
                      "Col1a1", "Pecam1")) +
  RotatedAxis() +
  ggtitle("Marker Expression Across Clusters")
dev.off()

# --- Assign labels ---
cat("\nAssigning cell type labels...\n")
cat("Cluster count:", length(levels(Idents(seurat_obj))), "\n")

new.cluster.ids <- c(
  "Olfactory Sensory Neurons",  # 0
  "Epithelial cells",           # 1
  "Macrophages",                # 2
  "Macrophages",                # 3
  "B cells",                    # 4
  "Olfactory Sensory Neurons",  # 5
  "Epithelial cells",           # 6
  "Neutrophils",                # 7
  "Macrophages",                # 8
  "T/NK cells",                 # 9
  "Fibroblasts",                # 10
  "Basal cells",                # 11
  "Epithelial cells",           # 12
  "Epithelial cells",           # 13
  "Epithelial cells",           # 14
  "Epithelial cells",           # 15
  "Epithelial cells",           # 16
  "Macrophages",                # 17
  "Sustentacular cells",        # 18
  "Olfactory Sensory Neurons",  # 19
  "Epithelial cells",           # 20
  "Epithelial cells",           # 21
  "Epithelial cells",           # 22
  "Dendritic cells",            # 23
  "B cells",                    # 24
  "Fibroblasts",                # 25
  "T cells",                    # 26
  "Macrophages",                # 27
  "Dendritic cells"             # 28
)

names(new.cluster.ids) <- levels(seurat_obj)
seurat_obj <- RenameIdents(seurat_obj, new.cluster.ids)
seurat_obj$cell_type <- Idents(seurat_obj)

pdf("plots/22_umap_annotated.pdf", width = 12, height = 8)
DimPlot(seurat_obj, reduction = "umap", label = TRUE, repel = TRUE) +
  ggtitle("Annotated UMAP") + NoLegend()
dev.off()

cat("Saving checkpoint...\n")
saveRDS(seurat_obj, "checkpoints/after_partB.rds")
cat("===== PART B COMPLETE =====\n")

RSCRIPT

echo "========== PART B DONE =========="
