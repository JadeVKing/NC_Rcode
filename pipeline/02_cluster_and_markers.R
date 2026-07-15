# =============================================================================
# GSE193816 scRNA-seq: Clustering and marker gene discovery
# =============================================================================

library(Seurat)
library(tidyverse)
library(qs)
library(RColorBrewer)
library(presto)

results_dir <- "results"

# ---- 1. Load clean data ----
seurat <- qread(file.path(results_dir, "seurat_clean_Ag&Pre.qs"))
DefaultAssay(seurat) <- "RNA"

# ---- 2. Elbow plot to assess PCs ----
ElbowPlot(seurat, ndims = 50)

# ---- 3. Cluster UMAP ----
clusters <- unique(seurat$seurat_clusters)
colors <- colorRampPalette(brewer.pal(10, "Set3"))(length(clusters))

p1 <- DimPlot(seurat, reduction = "umap", group.by = "seurat_clusters",
              label = TRUE, label.size = 4) +
  labs(title = "UMAP — Clusters") +
  scale_color_manual(values = setNames(colors, clusters)) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    text = element_text(size = 8, colour = "grey0"),
    axis.text = element_blank(),
    axis.title = element_blank(),
    legend.position = "right"
  )
ggsave(file.path(results_dir, "umap_clusters.png"), plot = p1, height = 4, width = 5, dpi = 300)

# Cluster composition
table(Idents(seurat))
table(seurat$sample, Idents(seurat))

# ---- 4. Variable features ----
VariableFeaturePlot(seurat)
top10 <- head(VariableFeatures(seurat), 10)

# ---- 5. Find all markers ----
markers <- FindAllMarkers(
  object = seurat,
  group.by = "seurat_clusters",
  assay = "RNA",
  only.pos = FALSE,
  logfc.threshold = 0.1,
  min.pct = 0.1,
  test.use = "wilcox"
)

write.csv(markers, file.path(results_dir, "global_markers.csv"), row.names = TRUE)

# Top 5 markers per cluster
top5 <- markers |>
  group_by(cluster) |>
  top_n(5, avg_log2FC)

# Heatmap
DoHeatmap(seurat,
  features = top5$gene,
  group.by = "seurat_clusters",
  slot = "scale.data",
  draw.lines = TRUE
)
ggsave(file.path(results_dir, "global_heatmap.png"), width = 8, height = 5)

# ---- 6. Save object ----
save(markers, seurat, file = file.path(results_dir, "seurat_ready.Rdata"))

message("Clustering and marker analysis complete. Markers saved to ", file.path(results_dir, "global_markers.csv"))
