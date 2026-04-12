#!/bin/bash
#SBATCH --job-name=scrnaseq_partD
#SBATCH --output=partD_%j.out
#SBATCH --error=partD_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4

module load r

echo "========== PART D: DE + Enrichment =========="

Rscript --no-save --no-restore - <<'RSCRIPT'

library(Seurat)
library(ggplot2)
library(DESeq2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)

dir.create("plots", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

cat("Loading checkpoint from Part B...\n")
seurat_obj <- readRDS("checkpoints/after_partB.rds")
cat("Loaded:", ncol(seurat_obj), "cells\n")

# --- Pseudobulk DE ---
cat("\nSetting up pseudobulk DE...\n")

seurat_obj_de <- subset(seurat_obj,
                        cells = colnames(seurat_obj)[
                          !is.na(seurat_obj$mouse_id) &
                          seurat_obj$mouse_id != ""])

macro_cells <- subset(seurat_obj_de, idents = "Macrophages")
macro_cells <- subset(macro_cells, subset = time %in% c("Naive", "D05"))
cat("Macrophages for DE:", ncol(macro_cells), "cells\n")

pseudo_macro <- AggregateExpression(macro_cells,
                                    assays        = "RNA",
                                    return.seurat = TRUE,
                                    group.by      = c("time", "mouse_id"))
Idents(pseudo_macro) <- "time"

de_results <- FindMarkers(pseudo_macro,
                          ident.1  = "D05",
                          ident.2  = "Naive",
                          test.use = "DESeq2")

cat("\nTop 20 DE genes:\n")
print(head(de_results, n = 20))

n_sig <- sum(de_results$p_val_adj < 0.05 & abs(de_results$avg_log2FC) > 1,
             na.rm = TRUE)
cat("\nSignificant genes (padj < 0.05, |log2FC| > 1):", n_sig, "\n")

write.csv(de_results, "results/de_Macrophages_D05_vs_Naive.csv")

# --- ORA ---
cat("\nRunning ORA...\n")

sig_genes <- rownames(de_results)[de_results$p_val_adj < 0.05 &
                                  abs(de_results$avg_log2FC) > 1]

if (length(sig_genes) > 0) {
  gene_ids <- bitr(sig_genes, fromType = "SYMBOL",
                   toType = "ENTREZID", OrgDb = org.Mm.eg.db)

  ora_results <- enrichGO(gene          = gene_ids$ENTREZID,
                          OrgDb         = org.Mm.eg.db,
                          ont           = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.05,
                          readable      = TRUE)

  pdf("plots/26_ora_dotplot.pdf", width = 10, height = 8)
  print(dotplot(ora_results, showCategory = 15) +
    ggtitle("GO Enrichment — Macrophages D05 vs Naive"))
  dev.off()

  pdf("plots/27_ora_barplot.pdf", width = 10, height = 8)
  print(barplot(ora_results, showCategory = 15))
  dev.off()

  write.csv(as.data.frame(ora_results), "results/ora_results.csv",
            row.names = FALSE)
  cat("ORA:", nrow(as.data.frame(ora_results)), "enriched terms\n")
} else {
  cat("No significant genes for ORA.\n")
}

# --- GSEA ---
cat("\nRunning GSEA...\n")

gene_list <- de_results$avg_log2FC
names(gene_list) <- rownames(de_results)
gene_list <- sort(gene_list, decreasing = TRUE)

gene_df <- bitr(names(gene_list), fromType = "SYMBOL",
                toType = "ENTREZID", OrgDb = org.Mm.eg.db)
gene_df <- gene_df[!duplicated(gene_df$SYMBOL), ]
gene_list_entrez <- gene_list[gene_df$SYMBOL]
names(gene_list_entrez) <- gene_df$ENTREZID
gene_list_entrez <- sort(gene_list_entrez, decreasing = TRUE)

gsea_results <- gseGO(geneList     = gene_list_entrez,
                      OrgDb        = org.Mm.eg.db,
                      ont          = "BP",
                      minGSSize    = 10,
                      maxGSSize    = 500,
                      pvalueCutoff = 0.05,
                      verbose      = FALSE)

if (nrow(as.data.frame(gsea_results)) > 0) {
  pdf("plots/28_gsea_dotplot.pdf", width = 10, height = 8)
  print(dotplot(gsea_results, showCategory = 15) +
    ggtitle("GSEA — Macrophages D05 vs Naive"))
  dev.off()

  pdf("plots/29_gsea_ridgeplot.pdf", width = 10, height = 8)
  print(ridgeplot(gsea_results, showCategory = 10))
  dev.off()

  n_sets <- min(3, nrow(as.data.frame(gsea_results)))
  pdf("plots/30_gsea_enrichment.pdf", width = 10, height = 8)
  print(gseaplot2(gsea_results, geneSetID = 1:n_sets))
  dev.off()

  write.csv(as.data.frame(gsea_results), "results/gsea_results.csv",
            row.names = FALSE)
  cat("GSEA:", nrow(as.data.frame(gsea_results)), "enriched gene sets\n")
} else {
  cat("No significant GSEA results.\n")
}

cat("===== PART D COMPLETE =====\n")

RSCRIPT

echo "========== PART D DONE =========="
