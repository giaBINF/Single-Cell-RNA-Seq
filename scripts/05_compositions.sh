#!/bin/bash
#SBATCH --job-name=scrnaseq_partE
#SBATCH --output=partE_%j.out
#SBATCH --error=partE_%j.err
#SBATCH --time=00:30:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4

module load r

echo "========== PART E: Composition + Extras =========="

Rscript --no-save --no-restore - <<'RSCRIPT'

library(Seurat)
library(dplyr)
library(ggplot2)

dir.create("plots", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

cat("Loading checkpoint from Part B...\n")
seurat_obj <- readRDS("checkpoints/after_partB.rds")

de_file <- "results/de_Macrophages_D05_vs_Naive.csv"
if (file.exists(de_file)) {
  de_results <- read.csv(de_file, row.names = 1)
  cat("Loaded DE results:", nrow(de_results), "genes\n")
} else {
  cat("WARNING: DE results not found. Run Part D first.\n")
  de_results <- NULL
}

# --- Composition by time ---
cat("\nGenerating composition plots...\n")

comp_data <- seurat_obj@meta.data %>%
  group_by(time, cell_type) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(time) %>%
  mutate(proportion = count / sum(count))

pdf("plots/31_composition_by_time.pdf", width = 10, height = 8)
ggplot(comp_data, aes(x = time, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_minimal() +
  labs(title = "Cell Type Composition Across Time Points",
       x = "Days Post Infection", y = "Proportion", fill = "Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# --- Composition by tissue and time ---
comp_data_tissue <- seurat_obj@meta.data %>%
  group_by(time, organ_custom, cell_type) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(time, organ_custom) %>%
  mutate(proportion = count / sum(count))

pdf("plots/32_composition_by_tissue_time.pdf", width = 14, height = 8)
ggplot(comp_data_tissue,
       aes(x = time, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~ organ_custom) +
  theme_minimal() +
  labs(title = "Cell Type Composition by Tissue and Time Point",
       x = "Days Post Infection", y = "Proportion", fill = "Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# --- Heatmap ---
if (!is.null(de_results)) {
  seurat_obj_de <- subset(seurat_obj,
                          cells = colnames(seurat_obj)[
                            !is.na(seurat_obj$mouse_id) &
                            seurat_obj$mouse_id != ""])

  pdf("plots/33_heatmap_de_genes.pdf", width = 12, height = 8)
  DoHeatmap(subset(seurat_obj_de, idents = "Macrophages"),
            features = head(rownames(de_results), 20)) +
    ggtitle("Top DE Genes in Macrophages (D05 vs Naive)")
  dev.off()
}

# --- Violin plots ---
pdf("plots/34_violin_isgs.pdf", width = 16, height = 8)
VlnPlot(seurat_obj,
        features = c("Isg15", "Ifit1", "Oasl2"),
        group.by = "cell_type",
        split.by = "time",
        pt.size  = 0)
dev.off()

pdf("plots/35_stacked_violin.pdf", width = 12, height = 8)
VlnPlot(seurat_obj,
        features = c("Isg15", "Ifit1", "Oasl2", "Cxcl10", "Ccl2"),
        stack    = TRUE,
        flip     = TRUE) +
  theme(legend.position = "right")
dev.off()

# --- Save final object ---
cat("Saving final annotated object...\n")
saveRDS(seurat_obj, "results/seurat_annotated.rds")

cat("\n===== PART E COMPLETE =====\n")
cat("===== FULL PIPELINE DONE =====\n")
cat("\nAll outputs:\n")
cat("  plots/     — 35 PDF figures\n")
cat("  results/   — DE tables, markers, GSEA, annotated .rds\n")

RSCRIPT

echo "========== PART E DONE =========="
