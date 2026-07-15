# =============================================================================
# GSE193816 scRNA-seq: Cell type annotation
# =============================================================================
# Uses manual annotation based on marker gene expression.
# Optional: SingleR / celldex for automated annotation (commented).
# =============================================================================

library(Seurat)
library(tidyverse)
library(qs)
library(RColorBrewer)
library(patchwork)

source("R/utils.R")

results_dir <- "results"

# ---- 1. Load data ----
load(file.path(results_dir, "seurat_ready.Rdata"))  # loads markers, seurat
# Alternative: seurat <- qread(file.path(results_dir, "seurat_clean_Ag&Pre.qs"))

# ---- 2. Marker genes for DotPlot preview ----
dotplot_genes <- c(
  "EPCAM", "KRT8", "KRT19",         # Epithelial
  "CD19", "CD79A", "MS4A1",         # B cells
  "MZB1", "JCHAIN", "GZMB",        # Plasma cells
  "CD3D", "CD3E", "CD4", "IL7R",    # T cells
  "CD8A", "CD8B",                   # CD8 T cells
  "NKG7", "GNLY",                   # NK cells
  "CD68", "LYZ", "CD14",            # Myeloid
  "CD1C", "CLEC4C",                  # cDC / pDC
  "TPSAB1", "CPA3", "KIT",          # Mast cells
  "FOXI1", "CFTR"                    # Ionocytes
)

# Preview annotation (cluster-level dotplot + UMAP)
result <- preview_annonation(seurat,
  dotplot_genes = dotplot_genes,
  output_prefix = "seurat",
  output_dir = results_dir
)

# ---- 3. Manual cell type mapping ----
# Adjust cluster numbers based on your data!
cell_type_map <- list(
  "Pulmonary ionocytes"       = c(16),
  "Airway epithelial cells"   = c(1, 2, 6, 7, 11, 15),
  "T Cells"                   = c(0, 3, 5),
  "B cells"                   = 8,
  "pDC"                       = 10,
  "cDC"                       = 4,
  "Neutrophils"               = 9,
  "Eosinophils"               = 14,
  "Macrophages"               = 12,
  "Mast cells"                = c(13)
)

seurat <- annotate_cell_types(seurat, cell_type_map, new_ident = "cell_type")

# Verify no missing clusters
table(seurat$cell_type, useNA = "ifany")
if (anyNA(seurat$cell_type)) {
  warning("Some clusters were not annotated! Check cell_type_map.")
}

# ---- 4. Reorder cell types ----
desired_order <- c(
  "B cells", "Airway epithelial cells", "cDC", "pDC",
  "Macrophages", "Pulmonary ionocytes", "Mast cells",
  "Neutrophils", "Eosinophils", "T Cells"
)
seurat$cell_type <- factor(seurat$cell_type, levels = desired_order)
Idents(seurat) <- "cell_type"

# ---- 5. Final annotated figures ----
annotated_final_plot(seurat,
  dotplot_genes = dotplot_genes,
  output_prefix = "seurat",
  output_dir = results_dir
)

# ---- 6. UMAP with cell type labels ----
umap_final <- CellDimPlot(seurat,
  label = TRUE, label_insitu = TRUE,
  reduction = "umap", group.by = "cell_type"
) +
  labs(title = "UMAP — Annotated Cell Types") +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        legend.position = "right") +
  guides(fill = "none")

ggsave(file.path(results_dir, "umap_annotated.png"),
  plot = umap_final, height = 5, width = 5.5, dpi = 300)

# ---- 7. DotPlot of marker genes ----
marker_genes <- c(
  "SCGB1A1", "KRT19", "AGR2",       # Epithelial
  "CD19", "MS4A1", "CD79B",          # B cells
  "CD1C", "CD1E", "CD1D",            # cDC
  "CLEC4C", "IL3RA", "SPIB",         # pDC
  "FOXI1", "CFTR", "SCNN1B",         # Ionocytes
  "CPA3", "KIT", "MS4A2",            # Mast cells
  "CXCL8", "S100A8", "S100A9",       # Neutrophils
  "CLC", "CCR3", "IL5RA",            # Eosinophils
  "CD3D", "CD3E", "NKG7"             # T cells
)

p_dot <- DotPlot(seurat,
  features = marker_genes,
  group.by = "cell_type",
  cols = c("lightgrey", "#E41A1C"),
  dot.scale = 5
) +
  labs(title = "Marker Gene Expression Across Cell Types") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1, colour = "black"),
    axis.text.y = element_text(size = 8, colour = "black"),
    plot.title = element_text(face = "bold", size = 10, hjust = 0.5)
  )

ggsave(file.path(results_dir, "dotplot_markers.png"),
  plot = p_dot, height = 5, width = 10, dpi = 300)

# ---- 8. Save annotated object ----
qsave(seurat, file.path(results_dir, "seurat_anno.qs"))

message("Cell type annotation complete. Saved to ", file.path(results_dir, "seurat_anno.qs"))
