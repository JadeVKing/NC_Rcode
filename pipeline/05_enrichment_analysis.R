# =============================================================================
# GSE193816 scRNA-seq: Pathway enrichment analysis (GSEA / GO / KEGG)
# =============================================================================
# Performs enrichment on differentially expressed genes between conditions.
# Note: Adjust OrgDb and species parameters for your dataset.
# =============================================================================

library(Seurat)
library(tidyverse)
library(qs)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)  # Change to org.Mm.eg.db for mouse

results_dir <- "results"

# ---- 1. Prepare DE gene list ----
# Load DE results or compute from seurat object
seurat <- qread(file.path(results_dir, "seurat_anno.qs"))
Idents(seurat) <- "condition"

de_results <- FindMarkers(seurat,
  ident.1 = "Asthma", ident.2 = "Control",
  only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25
)

# Filter significant genes
sig_genes <- de_results |>
  filter(p_val_adj < 0.05) |>
  arrange(desc(avg_log2FC))

# ---- 2. Prepare gene list for GSEA ----
# Convert gene symbols to ENTREZ IDs
gene_df <- bitr(rownames(sig_genes),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

gene_list <- sig_genes |>
  rownames_to_column("gene") |>
  inner_join(gene_df, by = c("gene" = "SYMBOL")) |>
  arrange(desc(avg_log2FC))

gsea_list <- gene_list$avg_log2FC
names(gsea_list) <- gene_list$ENTREZID
gsea_list <- sort(gsea_list, decreasing = TRUE)

# ---- 3. GSEA using GO terms ----
gse_go <- gseGO(
  geneList = gsea_list,
  OrgDb = org.Hs.eg.db,
  ont = "all",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

if (nrow(gse_go@result) > 0) {
  gse_go_res <- setReadable(gse_go, "org.Hs.eg.db", "ENTREZID")@result |>
    arrange(desc(enrichmentScore))
  write.csv(gse_go_res, file.path(results_dir, "GSEA_GO_results.csv"), row.names = TRUE)

  # Dotplot
  dotplot(gse_go, color = "p.adjust", showCategory = 10, split = ".sign") +
    facet_grid(. ~ .sign)
  ggsave(file.path(results_dir, "GSEA_GO_dotplot.pdf"), width = 8, height = 7)
}

# ---- 4. GSEA using KEGG ----
gse_kk <- gseKEGG(
  geneList = gsea_list,
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.9,
  verbose = FALSE
)

if (nrow(gse_kk@result) > 0) {
  gse_kk_res <- setReadable(gse_kk, "org.Hs.eg.db", "ENTREZID")@result |>
    arrange(desc(enrichmentScore))
  write.csv(gse_kk_res, file.path(results_dir, "GSEA_KEGG_results.csv"), row.names = TRUE)

  # GSEA plot for top pathway
  top_pathway <- gse_kk@result$ID[1]
  p <- gseaplot2(gse_kk, geneSetID = top_pathway,
    title = "", base_size = 10,
    pvalue_table = TRUE, rel_heights = c(1.5, 0.5, 1),
    ES_geom = "line"
  )
  ggsave(file.path(results_dir, paste0("GSEA_", top_pathway, ".pdf")),
    plot = p, width = 8, height = 6
  )
}

# ---- 5. Pathway activity scoring (GSVA / AddModuleScore) ----
if (exists("gse_kk_res") && nrow(gse_kk_res) > 0) {
  # Example: extract top pathway genes
  pathway_genes <- strsplit(gse_kk_res$core_enrichment[1], "/")[[1]]
  pathway_genes <- pathway_genes[pathway_genes %in% rownames(seurat)]

  if (length(pathway_genes) > 1) {
    seurat <- AddModuleScore(seurat,
      features = list(pathway_genes),
      name = "PathwayScore"
    )

    p_score <- VlnPlot(seurat, features = "PathwayScore1", group.by = "condition") +
      labs(title = gse_kk_res$Description[1])
    ggsave(file.path(results_dir, "pathway_activity_score.png"), plot = p_score,
           height = 3, width = 4, dpi = 300)
  }
}

message("Enrichment analysis complete. Results in ", results_dir, "/")
