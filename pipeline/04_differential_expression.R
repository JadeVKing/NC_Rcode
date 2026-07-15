# =============================================================================
# GSE193816 scRNA-seq: Differential expression and visualization
# =============================================================================

library(Seurat)
library(tidyverse)
library(qs)
library(ggsignif)
library(patchwork)

source("R/utils.R")

results_dir <- "results"

# ---- 1. Load annotated data ----
seurat <- qread(file.path(results_dir, "seurat_anno.qs"))

# ---- 2. Add condition metadata ----
# Derive condition from sample/Channel metadata
seurat$condition <- ifelse(grepl("ANA", seurat$Channel), "Control", "Asthma")
table(seurat$condition)

# ---- 3. DE analysis: Asthma vs Control ----
Idents(seurat) <- "condition"

global_markers_cond <- FindAllMarkers(seurat,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
write.csv(global_markers_cond, file.path(results_dir, "global_markers_by_condition.csv"),
          row.names = TRUE)

# ---- 4. Single gene exploration ----
# Example: PIEZO1 expression
markers_select <- FindMarkers(seurat,
  ident.1 = "Asthma", ident.2 = "Control",
  features = "PIEZO1"
)
print(markers_select)

# Violin plot
VlnPlot(seurat, features = c("PIEZO1"), group.by = "condition")

# Fetch expression data
expression_data <- FetchData(seurat,
  vars = c("PIEZO1", "IL6", "IL1B", "condition")
)

# Long format for plotting
long_data <- pivot_longer(expression_data,
  cols = -condition,
  names_to = "gene",
  values_to = "expression"
)

# Violin + boxplot for each gene
for (g in unique(long_data$gene)) {
  dat <- long_data |>
    filter(gene == g) |>
    mutate(condition = factor(condition, levels = c("Control", "Asthma")))

  p <- ggplot(dat, aes(x = condition, y = expression, fill = condition)) +
    geom_violin(trim = TRUE, alpha = 0.3, scale = "width", linewidth = 0.35) +
    geom_boxplot(alpha = 0.5, width = 0.5, outlier.size = 0.05) +
    geom_jitter(aes(color = condition), width = 0.2, size = 0.3, alpha = 0.8) +
    geom_signif(
      comparisons = list(c("Control", "Asthma")),
      map_signif_level = TRUE, test = "t.test",
      step_increase = 1, tip_length = 0, textsize = 3.5
    ) +
    labs(
      title = paste0(g, " Expression: Control vs Asthma"),
      y = bquote(italic(.(g)) ~ "Relative Expression")
    ) +
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line(linewidth = 0.35),
      axis.text = element_text(size = 8, colour = "grey0"),
      axis.title.x = element_blank(),
      plot.title = element_text(size = 8, face = "bold"),
      legend.position = "none"
    ) +
    scale_fill_manual(values = c("Control" = "royalblue", "Asthma" = "tomato")) +
    scale_color_manual(values = c("Control" = "royalblue", "Asthma" = "tomato"))

  ggsave(file.path(results_dir, paste0(g, "_violin.png")), plot = p, height = 3, width = 3, dpi = 300)
}

# ---- 5. Boxplot by cell type (grouped by condition) ----
genelist <- c("S100A9", "S100A8")
plots <- lapply(genelist, function(g) {
  boxplot_by_celltype(
    seurat_obj = seurat,
    gene = g,
    celltype_col = "cell_type",
    group_col = "condition",
    group_levels = c("Control", "Asthma"),
    color = c("Control" = "lightblue", "Asthma" = "tomato"),
    output_dir = results_dir,
    output_prefix = paste0(g, "_boxplot_by_celltype")
  )
})

message("DE analysis complete. Results in ", results_dir, "/")
