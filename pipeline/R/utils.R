# =============================================================================
# Utility functions for GSE193816 scRNA-seq analysis
# =============================================================================

# ---- Plot themes ----
base_theme <- theme_minimal() +
  theme(
    text = element_text(size = 8, face = "bold", hjust = 0.5),
    plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
    axis.title = element_blank(),
    axis.text = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none"
  )

text_theme <- theme(
  axis.text.x = element_text(angle = 45, hjust = 1),
  title = element_text(size = 8, colour = "grey0"),
  axis.line = element_line(linewidth = 0.35),
  axis.ticks = element_line(linewidth = 0.35, colour = "grey0"),
  text = element_text(size = 8, colour = "grey0"),
  axis.text = element_text(size = 8, colour = "grey0"),
  axis.title = element_text(size = 8, colour = "grey0")
)

# ---- Preview annotations with dotplot + UMAP ----
preview_annonation <- function(seurat_obj,
                               dotplot_genes,
                               output_prefix = "seurat",
                               output_dir = "results") {
  figure_dir <- file.path(output_dir, "preview_annotation")
  dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

  dotplot_colors <- c(low = "lightblue", high = "tomato")

  # DotPlot
  dotplot <- DotPlot(seurat_obj,
    features = dotplot_genes,
    group.by = "seurat_clusters",
    cols = dotplot_colors,
    dot.scale = 6
  ) + text_theme
  ggsave(file.path(figure_dir, paste0(output_prefix, "_dot_preview.png")),
    plot = dotplot, width = 12, height = 5, dpi = 300
  )

  # UMAP
  cell_types <- unique(seurat_obj$seurat_clusters)
  colors <- colorRampPalette(brewer.pal(10, "Set3"))(length(cell_types))
  umap_plot <- DimPlot(seurat_obj,
    reduction = "umap", group.by = "seurat_clusters",
    label = TRUE, label.size = 4, raster = TRUE
  ) +
    labs(title = "UMAP — Clusters") +
    scale_color_manual(values = setNames(colors, cell_types)) +
    base_theme
  ggsave(file.path(figure_dir, paste0(output_prefix, "_umap_preview.png")),
    plot = umap_plot, height = 6, width = 7, dpi = 300
  )

  invisible(list(seurat_obj = seurat_obj, dotplot = dotplot, umap_plot = umap_plot))
}

# ---- Cell type annotation ----
annotate_cell_types <- function(seurat_obj, cell_type_map, new_ident = "cell_type") {
  if (is.null(seurat_obj$seurat_clusters)) {
    seurat_obj$seurat_clusters <- Idents(seurat_obj)
  }
  cluster_to_celltype <- unlist(lapply(names(cell_type_map), function(celltype) {
    setNames(rep(celltype, length(cell_type_map[[celltype]])), cell_type_map[[celltype]])
  }))
  cluster_to_celltype <- setNames(as.character(cluster_to_celltype), names(cluster_to_celltype))
  annotation_vector <- cluster_to_celltype[as.character(seurat_obj$seurat_clusters)]
  names(annotation_vector) <- colnames(seurat_obj)
  seurat_obj <- AddMetaData(seurat_obj, metadata = annotation_vector, col.name = new_ident)
  message(sprintf("Added metadata column: %s", new_ident))
  seurat_obj
}

# ---- Final annotated plots ----
annotated_final_plot <- function(seurat_obj,
                                 dotplot_genes,
                                 output_prefix = "seurat",
                                 output_dir = "results") {
  figure_dir <- file.path(output_dir, "annotated_seurat")
  dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

  cell_types <- unique(seurat_obj$cell_type)
  colors <- colorRampPalette(brewer.pal(10, "Set3"))(length(cell_types))

  # UMAP
  umap_plot <- DimPlot(seurat_obj,
    reduction = "umap", group.by = "cell_type",
    label = TRUE, label.size = 4, raster = TRUE
  ) +
    labs(title = "UMAP — Annotated Cell Types") +
    scale_color_manual(values = setNames(colors, cell_types)) +
    base_theme

  ggsave(file.path(figure_dir, paste0(output_prefix, "_umap_anno.png")),
    plot = umap_plot, height = 6, width = 7, dpi = 300
  )

  # DotPlot
  dotplot <- DotPlot(seurat_obj,
    features = dotplot_genes,
    group.by = "cell_type",
    cols = c(low = "lightblue", high = "tomato"),
    dot.scale = 6
  ) + text_theme

  ggsave(file.path(figure_dir, paste0(output_prefix, "_dot_anno.png")),
    plot = dotplot, width = 10, height = 5, dpi = 300
  )

  # Pie chart
  df <- as.data.frame(table(seurat_obj$cell_type))
  pie_chart <- ggplot(df, aes(x = "", y = Freq, fill = Var1)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y") +
    scale_fill_manual(name = "Cell lineage", values = setNames(colors, unique(df$Var1))) +
    labs(title = "Cell Type Proportion") +
    theme_void() +
    theme(title = element_text(size = 8, colour = "grey0", hjust = 0.5))

  ggsave(file.path(figure_dir, paste0(output_prefix, "_pie.png")),
    plot = pie_chart, height = 5, width = 5, dpi = 300
  )

  # Cell proportion plots
  df_prop <- seurat_obj@meta.data |>
    group_by(cell_type, condition, sample) |>
    summarise(Count = n(), .groups = "drop") |>
    rename(Lineage = cell_type)

  p1 <- ggplot(df_prop, aes(
    y = fct_rev(fct_infreq(Lineage)),
    x = Count, fill = Lineage
  )) +
    geom_col() +
    labs(x = "Cell count", y = NULL) +
    scale_fill_manual(values = colors) +
    theme_bw(base_size = 10) +
    theme(axis.text = element_text(size = 8, color = "black"), legend.position = "none")

  p2 <- ggplot(df_prop, aes(x = Count, y = sample, fill = Lineage)) +
    geom_bar(position = "fill", stat = "identity") +
    labs(y = NULL, x = "Cell Lineage Proportion") +
    scale_fill_manual(name = "Cell lineage", values = colors) +
    theme_bw(base_size = 10) +
    theme(
      axis.text = element_text(size = 8, color = "black"),
      legend.key.size = unit(0.5, "cm")
    )

  coplot <- (p1 + p2) +
    plot_layout(guides = "collect") &
    plot_annotation(
      title = "Proportion of cell classification",
      theme = theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"))
    )

  ggsave(file.path(figure_dir, paste0(output_prefix, "_cell_proportion.png")),
    plot = coplot, height = 5, width = 8
  )

  invisible(list(
    seurat_obj = seurat_obj, dotplot = dotplot,
    umap_plot = umap_plot, coplot = coplot, pie_chart = pie_chart
  ))
}

# ---- Boxplot by cell type (grouped) ----
boxplot_by_celltype <- function(seurat_obj,
                                gene = "MERTK",
                                celltype_col = "cell_type",
                                group_col = "condition",
                                group_levels = c("Control", "Asthma"),
                                color = c("Control" = "lightblue", "Asthma" = "tomato"),
                                output_dir = "results",
                                output_prefix = "MERTK_boxplot") {
  figure_dir <- file.path(output_dir, "boxplots")
  dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

  if (!(gene %in% rownames(seurat_obj))) {
    stop(paste("Gene", gene, "not found."))
  }

  data <- FetchData(seurat_obj, vars = c(gene, celltype_col, group_col)) |>
    mutate(
      !!group_col := factor(.data[[group_col]], levels = group_levels),
      expression = .data[[gene]]
    )

  label_data <- data |>
    group_by(across(all_of(c(celltype_col, group_col)))) |>
    summarise(rate = mean(expression > 0) * 100, .groups = "drop")

  p <- ggplot(data, aes(x = .data[[celltype_col]], y = expression, fill = .data[[group_col]])) +
    geom_violin(trim = TRUE, alpha = 0.3, width = 0.8, scale = "width", linewidth = 0.35,
                position = position_dodge(0.8)) +
    geom_jitter(aes(color = .data[[group_col]]), size = 0.8, alpha = 0.6,
                position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.8)) +
    geom_text(data = label_data,
              aes(x = .data[[celltype_col]], y = -0.5,
                  label = paste0(round(rate, 1), "%"),
                  group = .data[[group_col]], color = .data[[group_col]]),
              position = position_dodge(width = 0.8),
              size = 3, show.legend = FALSE) +
    scale_fill_manual(values = color, guide = "none") +
    scale_color_manual(values = color) +
    labs(x = "Cell Subtype", y = bquote(italic(.(gene)) ~ "Relative Expression")) +
    theme_bw() +
    theme(
      legend.position = "top",
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 8, angle = 45, hjust = 1, colour = "black"),
      axis.text.y = element_text(size = 8, colour = "black"),
      panel.grid = element_blank(),
      panel.grid.major = element_line(linewidth = 0.35, colour = "grey95")
    )

  ggsave(file.path(figure_dir, paste0(output_prefix, ".png")), p, width = 6, height = 4, dpi = 300)
  ggsave(file.path(figure_dir, paste0(output_prefix, ".pdf")), p, width = 6, height = 4)

  p
}
