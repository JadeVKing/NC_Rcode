getwd()
rm(list = ls())
gc()
setwd("C:/wj/Ranalysis/sc_seq/GSE193816")

# setwd("C:/wj/Ranalysis/gse数据库/analysis")
# if (!requireNamespAEC("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("densmap")
# install.packages("DoubletFinder")
# remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
# library(devtools)
# install_github("mojaveazure/seurat-disk")
# remotes::install_github("immunogenomics/harmony")

# 单细胞质控功能模块 ----
source("C:/wj/R/Rscripts/seurat_toolbox.R")

if(T){
  register(SnowParam(workers = 4, progressbar = TRUE))
  setwd("C:/wj/Ranalysis/GEO/GSE193816/")

sceasy::convertFormat(
  obj = "GSE193816_all_data_raw_counts.h5ad",
  from = "anndata",
  to = "seurat",
  outFile = "scRNA.rds"
)

seurat <- readRDS("scRNA.rds")

DefaultAssay(seurat) <- "RNA"

colnames(seurat@meta.data)
head(seurat@meta.data)

# 先确认哪个列是 Ag / Pre / Dil
table(seurat$sample)      # 先看 sample 列
# table(seurat$condition) # 如果有 condition 也看
# table(seurat$group)     # 如果有 group 也看

# 假设 sample 列里确实是 Ag / Pre / Dil
seurat <- subset(seurat, subset = sample %in% c("Ag", "Pre"))
seurat$sample <- droplevels(seurat$sample)

table(seurat$sample)
  
  # Step 2: 质控
  # rownames(seurat) <- gsub("\\..*", "", rownames(seurat))
  
  seurat_obj <- basic_qc(seurat)
  
  seurat_obj <- filter_by_gene_umi(seurat_obj)

  seurat_obj <- subType(seurat_obj)
  
  # Step 4: 双细胞识别并移除
  df <- seurat_obj@meta.data
  table(seurat_obj$id)
  sce <- as.SingleCellExperiment(seurat_obj)
  sce <- scDblFinder(sce, samples = seurat_obj$Channel)
  seurat_obj$scDblFinder.class <- sce$scDblFinder.class
  rm(sce)
  gc()
  # 可视化双细胞结果（可选）
  CellDimPlot(seurat_obj, group.by = "scDblFinder.class", reduction = "umap")
  table(seurat_obj$scDblFinder.class)
  
  # Step 5: 移除 doublet 后重新聚类
  seurat_obj <- subset(seurat_obj, subset = scDblFinder.class == "singlet")
  # 切换到 RNA assay
  DefaultAssay(seurat_obj) <- "RNA"
  seurat_clean <- subType(seurat_obj,
                          res = 0.8,
                          sample_cells = 200000,
                          use_harmony = FALSE,
                          harmony_group = "sample_id")  # 重新运行
  
  seurat_clean <- FindNeighbors(seurat_clean, dims = 1:15)
  seurat_clean <- FindClusters(seurat_clean, resolution = 0.8)
  seurat_clean <- RunUMAP(seurat_clean, dims = 1:16, n.neighbors = 30)
  
  # Step 4: 保存干净数据----
  qsave(seurat_clean,"seurat_clean_Ag&Pre.qs")
  seurat_obj <- qread("seurat_clean_Ag&Pre.qs")
  # Step 5: 聚焦某类细胞群进行亚群分析（以 MNPs 为例）
  # MNPs <- subset(seurat, idents = "MNPs")
  # AEC_obj <- subType(MNPs)
  DimPlot(seurat_obj, group.by = "seurat_clusters", label = TRUE)
  # plot_markers(AEC_obj, marker_list)
}
gc()
#PCA 降维
ElbowPlot(seurat)
ElbowPlot(seurat, ndims = 50)
table(seurat_obj$seurat_clusters)

DimPlot(seurat_obj_Singlet, reduction = "umap", label = TRUE, group.by = "cell_set_clean") #如果样本差异很大说明有批次效应

seurat_clusters <- unique(seurat_obj$seurat_clusters)
colors <- colorRampPalette(brewer.pal(10, "Set3"))(length(seurat_clusters))
p1 <- DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE, label.size = 4) +
  labs(title = "UMAP OF seurat") +
  scale_color_manual(values = setNames(colors, seurat_clusters)) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    # title = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    text = element_text(size = 8, colour = "grey0"),  # 设置全局字体大小
    axis.text = element_blank(),
    axis.title = element_blank(),
    legend.position = "right") 
p1

# 保存
ggsave(paste0("global细胞类型", ".tiff"), plot = p1, height = 4, width = 5,dpi = 300)
ggsave(paste0("global细胞类", ".pdf"), plot = p1, height = 4, width = 5)
# 查看聚类结果
table(Idents(seurat_obj))
table(seurat_obj$sample, Idents(seurat_obj))
print(seurat_obj)

# 检查高变基因数量
length(VariableFeatures(seurat_obj))

# 可视化高变基因
VariableFeaturePlot(seurat_obj)

# 查看前 10 个变量基因
head(VariableFeatures(seurat_obj), 10)

# 标注前 10 个高变基因
top10 <- head(VariableFeatures(seurat_obj), 10)
# 绘制高变基因分布并标注前10个基因

# 绘制变量基因分布并标注
VariableFeaturePlot(seurat_obj) +
  geom_text_repel(data = data.frame(
    feature = top10,
    x = 1:length(top10),  # 可替换为实际的 PCA 或其他维度
    y = 1:length(top10)   # 可替换为实际的表达量或其他值
  ), aes(x = x, y = y, label = feature),
  max.overlaps = Inf)

# devtools::install_github('immunogenomics/presto')
# # 安装 presto 后，seurat_obj 会自动检测并使用 presto 的加速实现。
library(presto)
DefaultAssay(seurat_obj) <- "RNA"

# 为每个聚类寻找标记基因，验证是否存在生物学合理性：
markers <- FindAllMarkers(
  object = seurat_obj,
  group.by = "seurat_clusters",
  assay = "RNA",          # 根据预处理选择RNA或SCT
  only.pos = F,       # 保留上下调基因
  logfc.threshold = 0.1, # 过滤低差异基因
  min.pct = 0.1,         # 基因在至少10%的细胞中表达
  test.use = "wilcox"     # 默认Wilcoxon秩和检验
)
head(markers)
write.csv(markers,file = "global_markers.csv",row.names = T)
colnames(seurat_obj@meta.data)
# 统计聚类分组中的细胞数量
table(seurat_obj$seurat_obj_clusters)
# 统计样本分组中的细胞数量
table(seurat_obj$sample)

# 统计自定义分组中的细胞数量
table(seurat_obj$percent.ribo)
# 数据备份----
save(markers,seurat_obj,file = "seurat_ready.Rdata")
load("seurat_ready.Rdata")

seurat_obj <- readRDS("seurat_Ag&Pre.rds")
Layers(seurat_obj)
DefaultAssay(seurat_obj)  # 确认当前活跃 assay
gc()

GetAssayData(seurat_obj, layer  = "data")[1:5, 1:5]       # 查看标准化后的数据
GetAssayData(seurat_obj, layer  = "scale.data")[1:5, 1:5] # 查看缩放后的数据
top5 <- markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)
# 使用标准化后的数据（data 层）绘制基因表达热图
DoHeatmap(seurat_obj, 
          label = T,
          features = top5$gene, 
          group.by = "cell_set_clean",
          slot = "scale.data",
          draw.lines = T)
# Save plot
ggsave(
  filename = "GLOBALDoHeatmap图.pdf",
  width = 8,  # Slightly wider to accommodate legend
  height = 5   # Adjusted aspect ratio
)
ggsave(
  filename = "GLOBALDoHeatmap.png",
  width = 8,  # Slightly wider to accommodate legend
  height = 5   # Adjusted aspect ratio
)

# 差异表达分析----
Idents(seurat_obj) <- "seurat_clusters"
global_markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(global_markers, file = "global_markers1013.csv", row.names = T)

table(seurat_obj$phenotype)
table(seurat_obj$id)
table(seurat_obj$sample)
table(seurat_obj$Channel)
table(seurat_obj$condition)

# 添加分组信息----
seurat_obj$condition <- ifelse(grepl("ANA", seurat_obj$Channel), "Control", "Asthma")
table(seurat_obj$condition)
# 设置分组身份，差异分析将默认此分组
Idents(seurat_obj) <- "condition"


cell_types <- unique(seurat_obj$cell_set_clean)
colors <- colorRampPalette(brewer.pal(10, "Set3"))(length(cell_types))
# colors <- c("0"="#8DD3C7", "1"="#ECEBBE", "2"="#E0989E", "3"="#91AAC5", 
#             "4"="#F2BA63", "5"="#D2D69E", "6"="#E3D5DC", "7"="#BC80BD")
# Create cell count plot (left panel)
# 条形图----
if(T){
  # cell proportion
  df <- seurat_obj@meta.data
  df_prop <- df %>%
    dplyr::group_by(cell_set_clean,condition, Channel, sample) %>%
    dplyr::summarise(Count = dplyr::n(), .groups = "drop") %>% 
    as.data.frame() %>% 
    setNames(c("Lineage", "condition", "sample","state","Count"))
  
  df_prop <- df_prop %>%
    mutate(
      condition = factor(condition, levels = c("Control", "Asthma")),
      sample = factor(sample, levels = unique(sample[order(condition)]))  # 按condition排序
    )
  df_prop <- arrange(df_prop, condition)
  
  p1 <- df_prop %>%
    ggplot(aes(y = forcats::fct_rev(forcats::fct_infreq(Lineage)), 
               x = Count,  # 使用预计算的 Count 列
               fill = Lineage)) +
    geom_col() +  # 用 geom_col 直接绘制列高
    labs(x = 'Cell count', y = NULL) +
    scale_fill_manual(values = colors) +
    theme_bw(base_size = 10) +
    theme(
      axis.text = element_text(size = 8, color = 'black'),
      legend.position = "none"
    )
  
# Create frequency plot (right panel)
p2 <- ggplot(df_prop, 
             aes(x = Count,
                 y = sample,
                 fill = Lineage)) +
  geom_bar(position = "fill", stat = "identity") +
  # coord_flip() +
  labs(y = NULL, x = 'Cell Lineage Proportion') +
  # scale_x_discrete(labels = sample_labels) +
  scale_fill_manual(
    name = "Cell lineage",
    values = colors
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text = element_text(size = 8, color = 'black'),
    legend.key.size = unit(0.5, "cm")  # Adjust legend key size
  ) 
library(patchwork)
# Combine plots
coplot <- (p1 + p2) + 
  plot_layout(guides = "collect") & 
  plot_annotation(
    title = "Proportion of cell classification",
    theme = theme(
      plot.title = element_text(
        hjust = 0.5, 
        size = 10, 
        face = 'bold'
      )
    )
  )
}
coplot


# Save plot
ggsave(
  filename = "GLOBAL分组细胞比例图.pdf",
  plot = coplot,
  width = 8,  # Slightly wider to accommodate legend
  height = 5   # Adjusted aspect ratio
)
ggsave(
  filename = "GLOBAL分组细胞比例图.png",
  plot = coplot,
  width = 8,  # Slightly wider to accommodate legend
  height = 5   # Adjusted aspect ratio
)

# 差异表达分析
markers_select <- FindMarkers(seurat_obj, ident.1 = "Asthma", ident.2 = "Control", features = "PIEZO1")
markers_select <- FindMarkers(seurat_obj, ident.1 = "AA", ident.2 = "ANA")
print(markers_select)
VlnPlot(seurat_obj_Ag, features = c("PIEZO1"), group.by = "phenotype") #初步查看

# 默认从 `data` slot 提取标准化数据
expression_data <- FetchData(seurat_obj, vars = c("PIEZO1"))

# 从 `counts` slot 提取原始计数数据
expression_data_raw <- FetchData(seurat_obj, vars = c("PIEZO1"), slot = "counts")

# 从 `scale.data` slot 提取 Z-score 标准化数据
expression_data_zscore <- FetchData(seurat_obj, vars = c("PIEZO1"), slot = "scale.data")


if(T){
  expression_data <- as.data.frame(FetchData(seurat_obj_Ag, vars = c("PIEZO1","IL6","IL1B","phenotype")))
  
  library(tidyr)
  long_predata <- pivot_longer(
    data = expression_data,
    cols = !c(phenotype),      # 匹配以 "Gender_" 开头的列
    names_to = "name",             # 新的变量列名称
    values_to = "value"            # 新的值列名称
  )
  table(long_predata$name)
  # # 显示转换后的长数据格式的数据框
  print(long_predata)
  # 检查数据
  head(expression_data)
  for (name in unique(long_predata$name)){
    expression_data <- long_predata[grep(paste0("^",name,"$"),long_predata$name, ignore.case = TRUE), ] 
    
    # 小提琴图
    expression_data <- expression_data %>% mutate(group = factor(phenotype, levels = c("ANA", "AA")))
    expression_data <- arrange(expression_data, group)
    # 移除表达为0的样本
    # expression_data <- expression_data[expression_data$value > 0, ]
    
    # 绘制箱线图
    p <- ggplot(expression_data, aes(x = group, y = value, fill = group)) +
      geom_violin(trim = T, alpha = 0.3, scale = "width", linewidth = 0.35) + # 增加小提琴图，透明度调低
      geom_boxplot(aes(fill = group),
                   alpha = 0.5,
                   linewidth = 0.5,
                   notch = FALSE,
                   width = 0.5,
                   outlier.size = 0.05,
                   outlier.color = "black"
      ) +
      geom_jitter(mapping = aes(colour = group),
                  shape = 16,
                  position = position_jitter(0.2),
                  size = 0.3,
                  alpha = 0.8
      ) +
      geom_signif(
        comparisons = list(c("ANA", "AA")), 
        map_signif_level = T, 
        test = "t.test", # 或 "wilcox.test" 视具体情况
        step_increase = 1,
        tip_length = 0,
        textsize = 3.5,
        size = 0.35
      ) +
      stat_summary(fun = mean, geom = "point", color = "blue", size = 2) +
      stat_summary(fun = median, geom = "point", color = "red", size = 2) +
      labs(
        title = paste0(name, " Expression in Healthy vs Disease"), 
        x = "Condition", 
        y = bquote(italic(.(name)) ~ " Relative Expression")
      ) +
      theme(
        # panel.grid.major = element_line(color = "grey90", size = 0.35),
        # panel.border = element_rect(color = "black", fill = NA, size = 0.35),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(linewidth = 0.35),
        axis.ticks = element_line(linewidth = 0.35, colour = "grey0"),
        text = element_text(size = 8, colour = "grey0"),
        axis.text.x = element_text(size = 8,  colour = "grey0"),
        axis.text.y = element_text(size = 8, colour = "grey0"),
        axis.title.y = element_text(size = 8, colour = "grey0"),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 8, fAEC = "bold"),
        legend.position = "none"
      ) +
      scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
      # scale_y_log10(expand = expansion(mult = c(0, 0.1))) +
      scale_fill_manual(values = c("ANA" = "royalblue", "AA" = "tomato")) +
      scale_colour_manual(values = c("ANA" = "royalblue", "AA" = "tomato"))
    # 显示图
    print(p)
    # 保存图片,多个柱子，宽度加0.5
    ggsave(paste0(name, ".tiff"), plot = p, height = 3, width = 3,dpi = 300)
    ggsave(paste0(name, ".pdf"), plot = p, height = 3, width = 3,dpi = 300)
    # 
    # result <- t.test(value~phenotype, expression_data, var.equal = T)#方差齐
    # print(result)
  }
}

# 细胞定义----

preview_annonation <- function(seurat_obj,
                               dotplot_genes,
                               output_prefix = "seurat",
                               save_rds = F,
                               output_dir = ".") {
  
  figure_dir <- file.path(output_dir, "preview_annonation")
  dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 设置颜色
  dotplot_colors <- c(low = "lightblue", high = "tomato")
  
  
  # DotPlot 绘图----
  dotplot <- DotPlot(seurat_obj,
                     features = dotplot_genes,
                     group.by = "seurat_clusters",
                     cols = dotplot_colors,
                     dot.scale = 6) +
    text_theme
  ggsave(file.path(figure_dir, paste0(output_prefix, "_dot_anno.png")),
         plot = dotplot, width = 12, height = 5, dpi = 300)
  ggsave(file.path(figure_dir, paste0(output_prefix, "_dot_anno.pdf")),
         plot = dotplot, width = 12, height = 5, dpi = 300)
  
  
  # 分配注释----
  cell_types <- unique(seurat_obj$seurat_clusters)
  colors <- colorRampPalette(brewer.pal(10, "Set3"))(length(cell_types))
  
  # UMAP 图
  umap_plot <- DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters",
                       label = TRUE, label.size = 4, raster = TRUE) +
    labs(title = "UMAP with Annotated Cell Types") +
    scale_color_manual(values = setNames(colors, cell_types)) +
    base_theme
  
  ggsave(file.path(figure_dir, paste0(output_prefix, "_umap.png")),
         plot = umap_plot, height = 6, width = 7, dpi = 300)
  ggsave(file.path(figure_dir, paste0(output_prefix, "_umap.pdf")),
         plot = umap_plot, height = 6, width = 7)
  
  
  return(list(seurat_obj = seurat_obj,
              dotplot = dotplot,
              umap_plot = umap_plot))
}

annotate_cell_types <- function(seurat_obj, cell_type_map, new_ident = "cell_type") {
  
  # 确保 seurat_clusters 已定义
  if (is.null(seurat_obj$seurat_clusters)) {
    seurat_obj$seurat_clusters <- Idents(seurat_obj)
  }
  
  # 构建映射向量（cluster ID -> cell type）
  cluster_to_celltype <- unlist(lapply(names(cell_type_map), function(celltype) {
    setNames(rep(celltype, length(cell_type_map[[celltype]])), cell_type_map[[celltype]])
  }))
  cluster_to_celltype <- setNames(as.character(cluster_to_celltype), names(cluster_to_celltype))
  
  # 创建一个新的 metadata 列，按 cluster 映射 cell type
  annotation_vector <- cluster_to_celltype[as.character(seurat_obj$seurat_clusters)]
  names(annotation_vector) <- colnames(seurat_obj)
  
  # 添加到 Seurat 对象中
  seurat_obj <- AddMetaData(seurat_obj, metadata = annotation_vector, col.name = new_ident)
  
  # 保留原有的 Idents，不进行更改
  if (new_ident %in% colnames(seurat_obj@meta.data)) {
    message(sprintf("Added new metadata column: %s", new_ident))
  }
  
  return(seurat_obj)
}


annotated_final_plot <- function(seurat_obj,
                                 dotplot_genes,
                                 output_prefix = "seurat",
                                 save_qs = F,
                                 output_dir = ".") {
  
  figure_dir <- file.path(output_dir, "annotated_seurat")
  dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
  
  cell_types <- unique(seurat_obj$cell_type)
  colors <- colorRampPalette(brewer.pal(10, "Set3"))(length(cell_types))
  
  # UMAP 图
  umap_plot <- DimPlot(seurat_obj, reduction = "umap", group.by = "cell_type",
                       label = TRUE, label.size = 4, raster = TRUE) +
    labs(title = "UMAP with Annotated Cell Types") +
    scale_color_manual(values = setNames(colors, cell_types)) +
    base_theme
  
  ggsave(file.path(figure_dir, paste0(output_prefix, "_umap.png")),
         plot = umap_plot, height = 6, width = 7, dpi = 300)
  ggsave(file.path(figure_dir, paste0(output_prefix, "_umap.pdf")),
         plot = umap_plot, height = 6, width = 7)
  
  # 设置颜色
  dotplot_colors <- c(low = "lightblue", high = "tomato")
  # scale_color_gradient(low = dotplot_colors["low"], high = dotplot_colors["high"])
  text_theme <- theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    title = element_text(size = 8, colour = "grey0"),
    axis.line = element_line(linewidth = 0.35),
    axis.ticks = element_line(linewidth = 0.35, colour = "grey0"),
    text = element_text(size = 8, colour = "grey0"),
    axis.text = element_text(size = 8, colour = "grey0"),
    axis.title = element_text(size = 8, colour = "grey0"))
  
  # DotPlot 绘图----
  dotplot <- DotPlot(seurat_obj,
                     features = dotplot_genes,
                     group.by = "cell_type",
                     cols = dotplot_colors,
                     dot.scale = 6) +
    text_theme
  ggsave(file.path(figure_dir, paste0(output_prefix, "_dot_anno.png")),
         plot = dotplot, width = 10, height = 5, dpi = 300)
  ggsave(file.path(figure_dir, paste0(output_prefix, "_dot_anno.pdf")),
         plot = dotplot, width = 10, height = 5, dpi = 300)
  
  # 饼图----
  df <- as.data.frame(table(seurat_obj$cell_type))
  pie_chart <- ggplot(df, aes(x = "", y = Freq, fill = Var1)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y") +
    scale_fill_manual(name = "Cell lineage", values = setNames(colors, unique(df$Var1))) +
    labs(title = "Cell Type Proportion") +
    theme_void() +
    theme(
      title = element_text(size = 8, colour = "grey0", hjust = 0.5),
      text = element_text(size = 8, colour = "grey0"))
  
  ggsave(file.path(figure_dir, paste0(output_prefix, "_pie.tiff")),
         plot = pie_chart, height = 5, width = 5, dpi = 300)
  ggsave(file.path(figure_dir, paste0(output_prefix, "_pie.pdf")),
         plot = pie_chart, height = 5, width = 5)
  
  # 条形图----
  
  get_df_prop <- function(seurat_obj) {
    df <- seurat_obj@meta.data
    df_prop <- df %>%
      dplyr::group_by(cell_set_clean, phenotype, Channel, sample) %>%
      dplyr::summarise(Count = dplyr::n(), .groups = "drop") %>% 
      as.data.frame() %>%
      setNames(c("Lineage", "condition", "sample", "state", "Count")) %>%
      mutate(
        condition = factor(condition, levels = c("AA", "ANA")),
        sample = factor(sample, levels = unique(sample[order(condition)]))
      ) %>%
      arrange(condition)
    return(df_prop)
  }
  
  df_prop <- get_df_prop(seurat_obj)
  
  p1 <- df_prop %>%
    ggplot(aes(y = forcats::fct_rev(forcats::fct_infreq(Lineage)), 
               x = Count,  # 使用预计算的 Count 列
               fill = Lineage)) +
    geom_col() +  # 用 geom_col 直接绘制列高
    labs(x = 'Cell count', y = NULL) +
    scale_fill_manual(values = colors) +
    theme_bw(base_size = 10) +
    theme(
      axis.text = element_text(size = 8, color = 'black'),
      legend.position = "none"
    )
  
  # Create frequency plot (right panel)
  p2 <- ggplot(df_prop, 
               aes(x = Count,
                   y = sample,
                   fill = Lineage)) +
    geom_bar(position = "fill", stat = "identity") +
    # coord_flip() +
    labs(y = NULL, x = 'Cell Lineage Proportion') +
    # scale_x_discrete(labels = sample_labels) +
    scale_fill_manual(
      name = "Cell lineage",
      values = colors
    ) +
    theme_bw(base_size = 10) +
    theme(
      axis.text = element_text(size = 8, color = 'black'),
      legend.key.size = unit(0.5, "cm")  # Adjust legend key size
    ) 
  library(patchwork)
  # Combine plots
  coplot1 <- (p1 + p2) + 
    plot_layout(guides = "collect") & 
    plot_annotation(
      title = "Proportion of cell classification",
      theme = theme(
        plot.title = element_text(
          hjust = 0.5, 
          size = 10, 
          face = 'bold'
        )
      )
    )
  
  coplot1
  
  ggsave(file.path(figure_dir, paste0(output_prefix, "细胞比例堆叠图.png")),
         plot = coplot1, height = 5, width = 8)
  ggsave(file.path(figure_dir, paste0(output_prefix, "细胞比例堆叠图.pdf")),
         plot = coplot1, height = 5, width = 8)
  
  # 条形图----
  p3 <- df_prop %>%
    ggplot(aes(y = forcats::fct_rev(forcats::fct_infreq(Lineage)), 
               x = Count,  # 使用预计算的 Count 列
               fill = Lineage)) +
    geom_col() +  # 用 geom_col 直接绘制列高
    labs(x = 'Cell count', y = NULL) +
    scale_fill_manual(values = colors) +
    theme_bw(base_size = 10) +
    theme(
      axis.text = element_text(size = 8, color = 'black'),
      legend.position = "none"
    )
  
  # Create frequency plot (right panel)
  p4 <- ggplot(df_prop, 
               aes(x = Count,
                   y = condition,
                   fill = Lineage)) +
    geom_bar(position = "fill", stat = "identity") +
    # coord_flip() +
    labs(y = NULL, x = 'Cell Lineage Proportion') +
    # scale_x_discrete(labels = sample_labels) +
    scale_fill_manual(
      name = "Cell lineage",
      values = colors
    ) +
    theme_bw(base_size = 10) +
    theme(
      axis.text = element_text(size = 8, color = 'black'),
      legend.key.size = unit(0.5, "cm")  # Adjust legend key size
    ) 
  library(patchwork)
  # Combine plots
  coplot2 <- (p3 + p4) + 
    plot_layout(guides = "collect") & 
    plot_annotation(
      title = "Proportion of cell classification",
      theme = theme(
        plot.title = element_text(
          hjust = 0.5, 
          size = 10, 
          face = 'bold'
        )
      )
    )
  
  coplot2
  
  ggsave(file.path(figure_dir, paste0(output_prefix, "分组细胞比例堆叠图.png")),
         plot = coplot2, height = 5, width = 8)
  ggsave(file.path(figure_dir, paste0(output_prefix, "分组细胞比例堆叠图.pdf")),
         plot = coplot2, height = 5, width = 8)
  
  # 保存 qs----
  if (save_qs) {
    qsave(seurat_obj, file = paste0(output_prefix, "_anno.qs"))
  }
  return(list(seurat_obj = seurat_obj,
              dotplot = dotplot,
              umap_plot = umap_plot,
              coplot1 = coplot1,
              coplot2 = coplot2,
              pie_chart = pie_chart))
}

dotplot_genes <- c(
  "LYZ","S100A12", "APOBEC3A","CSF3R", # 炎症单核细胞（Classical monocytes）
  "MARCO", "APOE","MRC1","CD163","MSR1", #M2
  "CD80", "CD86","NOS2","IL1B","IL6", #M1
  "UBE2C", "RRM2", "MKI67", "CDC20", #增殖巨噬细胞细胞群
  "FCER1A","HLA-DQB2", "CLEC10A", "CD1C","CD1E", #DC2
  "TLR9", "CLEC4C", "GZMB", "VASH2", #pDC
  # 成纤维细胞
  "COL1A1", "COL1A2", "COL3A1", "LUM", "DCN",
  # 肌成纤维
  "COL6A1","MYLK","CALD1","PDGFRB",
  # Basal
  "TP63", "KRT5", "KRT17", "KRT15", "DSC3",
  # # Cycling
  # "MKI67", "TYMS", "CENPF", "TOP2A", "RRM2",
  # Club
  "SCGB1A1", "SCGB3A1", "SCGB3A2", "SFTA1P", "RNASE1",
  # Ciliated
  "FOXJ1", "PIFO", "DNAH5", "TUBA1A", "RSPH4A",
  # AT2
  "SFTPC", "SFTPB", "SFTPD", "NAPSA", "SFTPA1",
  # AT1
  "AGER", "RTKN2", "SPOCK2", "CAV1","CAV2",
  # Goblet
  "ZG16B", "MUC5B", "VMO1", "TFF3", "BPIFB1",
  # 内皮细胞
  "PECAM1", "VWF", "CDH5", "CLDN5","SPARCL1",
  # 肥大细胞
  "TPSAB1", "TPSB2", "TPSD1", "CPA3", "KIT",
  # T细胞
  "CD3D", "CD3E", "CD2", "TRAC", "IL7R", "CD4", "CD8A",
  # B细胞
  "CD19", "CD79A", "MS4A1", "CD79B", "BANK1",
  # Plasma
  "MZB1", "XBP1", "SDC1", "IGHG1", "PRDM1",
  # # NK细胞
  "NKG7", "GNLY", "PRF1", "KLRB1", "KLRD1",
  # # 中性粒细胞
  "S100A8", "S100A9", "MPO", "FCGR3B","MALT1","MIF",
  # smooth muscle cells
  "ACTA2","MYH11","TAGLN","TPM1","MYL9",
  # 间皮细胞
  "WT1", "MSLN", "CALB2", "UPK3B", "PDPN" 
)

# 查看元数据
df <- seurat@meta.data

# 设定主题----
# 1.适合umap，无坐标
base_theme <- theme_minimal() +  #  theme_classic()
  theme(
    text = element_text(size = 8, face = "bold", hjust = 0.5),
    plot.title = element_text(size = 8, face = "bold", hjust = 0.5),  # 标题居中 & 加粗
    axis.title = element_blank(),  # 隐藏坐标轴标题
    axis.text = element_blank(),
    panel.grid = element_blank(),  # 去除网格线
    legend.position = "none"      # 图例放右侧（可改 "bottom", "top"）
  )
# 2.适合dotplot，有坐标
text_theme <- theme(
  axis.text.x = element_text(angle = 45, hjust = 1),
  title = element_text(size = 8, colour = "grey0"),
  axis.line = element_line(linewidth = 0.35),
  axis.ticks = element_line(linewidth = 0.35, colour = "grey0"),
  text = element_text(size = 8, colour = "grey0"),
  axis.text = element_text(size = 8, colour = "grey0"),
  axis.title = element_text(size = 8, colour = "grey0"))

# 注释----
seurat_obj$cell_type <- seurat_obj$cell_set_clean
result <- preview_annonation(seurat_obj,
                             dotplot_genes = dotplot_genes,
                             output_prefix = "seurat",
                             save_rds = F,
                             output_dir = "results")
result
rm("result")
gc()
# 重新定义
cell_type_map <- list(
  "Pulmonary ionocytes" = c(16),
  "Airway epithelial cells" = c(1,2,6,7,11,15),
  "T Cells" = c(0,3,5),
  "B cells" = 8,
  "pDC" = 10,
  "cDC" = 4,
  "Neutrophils" = 9,
  "Eosinophils" = 14,
  "Macrophages" = 12,
  "Mast cells" = c(13)
)

seurat_obj <- annotate_cell_types(seurat_obj, cell_type_map)

# 查看是否定义漏了----
# 查看 meta.data 是否有 cell_type 列
if(!"cell_type" %in% colnames(seurat_obj@meta.data)){
  stop("seurat_obj@meta.data 中没有 cell_type 这一列！")
}

# 查看 cell_type 列的值分布
table(seurat_obj$cell_type, useNA = "ifany")

# 出图----
# SCP绘制UMAP----
umap_plot <- CellDimPlot(seurat_obj, 
                         label = T,
                         label_insitu = T,
                         reduction = "umap", 
                         group.by = "cell_type") +
  labs(title = "UMAP of Seurat Object") +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),  # 标题加粗、居中，字号稍大
    legend.position = "right"
  ) +
  guides(fill = "none") 
umap_plot
ggsave("umap1017.png",
       plot = umap_plot, height = 5, width = 5.5, dpi = 300)
ggsave("umap1017.pdf",
       plot = umap_plot, height = 5, width = 5.5, dpi = 300)

annotated_final_plot(seurat_obj,
                     dotplot_genes = dotplot_genes,
                     output_prefix = "seurat",
                     save_qs = F,
                     output_dir = ".")
cell_type_map <- list(
  "Pulmonary ionocytes" = c(16),
  "Airway epithelial cells" = c(1,2,6,7,11,15),
  "T Cells" = c(0,3,5),
  "B cells" = 8,
  "pDC" = 10,
  "cDC" = 4,
  "Neutrophils" = 9,
  "Eosinophils" = 14,
  "Macrophages" = 12,
  "Mast cells" = c(13)
)

desired_order <- c("B cells", "Airway epithelial cells", "cDC", "pDC", "Macrophages", "Pulmonary ionocytes","Mast cells", "Neutrophils", "Eosinophils","T Cells")
seurat_obj$cell_type <- factor(seurat_obj$cell_type, levels = desired_order)
Idents(seurat_obj) <- "cell_type"

# DotPlot 绘图----
library(Seurat)
library(ggplot2)

dotplot_genes <- c(
  # 上皮细胞
  "EPCAM", "KRT8", "KRT19",
  # B细胞
  "CD19", "CD79A", "MS4A1", 
  # 浆细胞
  "MZB1", "JCHAIN","GZMB", 
  # T细胞
  "CD3D", "CD3E", 
  "CD4", "IL7R",     # CD4+ T
  "CD8A", "CD8B",    # CD8+ T
  # 肥大细胞
  "TPSAB1", "CPA3", "KIT",
  # 单核-巨噬细胞系统
  "CD68", "LYZ", "CD14",
  # NK细胞
  "NKG7", "GNLY", "KLRB1"
)
dotplot_genes <- c(
  # # AT1
  # "AGER", "RTKN2", "CAV1",
  # # AT2
  # "SFTPC", "SFTPB", "SFTPD",
  # B细胞
  "CD19", "MS4A1", "CD79B",
  # Club
  "SCGB1A1", "KRT19", "AGR2",
  # Basal
  # "KRT5", "KRT17", "KRT15",
  # # Cycling
  # "MKI67", "TYMS", "CENPF", "TOP2A", "RRM2",
  # Ciliated,Secretory cells
  # "PIFO", "DNAH5", "TUBA1A",
  # Goblet
  # "ZG16B", "MUC5B", "VMO1", "TFF3", "BPIFB1",
  #cDC
  "CD1C","CD1E","CD1D",
  # pDC
  "CLEC4C","IL3RA","SPIB",
  # 成纤维细胞
  # "COL1A1", "COL1A2", "COL3A1",
  #M2
  "MARCO","MSR1", "LYZ", 
  # pulmonary ionocyte
  "FOXI1","CFTR","SCNN1B",
  # 肥大细胞
  "CPA3", "KIT","MS4A2",
  # 间皮细胞
  # "WT1", "UPK3B", "PDPN",
  # 炎症单核细胞（Classical monocytes）
  # "CSF1R", "AZU1", "CD14", 
  # # 中性粒细胞
  "CXCL8", "S100A8", "S100A9",
  # 嗜酸性粒细胞
  "CLC","CCR3","IL5RA",
  # smooth muscle cells
  # "ACTA2","MYH11","TAGLN",
  # T细胞
  "CD3D", "CD3E", "NKG7"
)
# 绘制 DotPlot
p <- DotPlot(seurat_obj,
             features = dotplot_genes,
             group.by = "seurat_clusters",
             cols = c("lightgrey", "#E41A1C"), # 用更柔和的红色
             dot.scale = 5) + 
  labs(
    x = NULL,
    y = NULL,
    title = "Marker Gene Expression Across Cell Populations"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major = element_line(linewidth = 0.3, colour = "grey90"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1, vjust = 1, colour = "black"),
    axis.text.y = element_text(size = 8, colour = "black"),
    axis.title = element_text(face = "bold", size = 9),
    legend.title = element_text(face = "bold", size = 9),
    legend.text = element_text(size = 8),
    plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
    aspect.ratio = 0.5  # 让点图更扁一些，节省横向空间
  )

# 如果希望更明显区分不同细胞类型，可以加 facet 或手动调整顺序
# seurat_obj$cell_type <- factor(seurat_obj$cell_type, levels = c("Epithelial", "Macrophage", ...))
p

ggsave("1013cellmarker_dotplot.png",
       plot = p, height = 8, width = 16, dpi = 300)
ggsave("1013cellmarker_dotplot.pdf",
       plot = p, height = 5, width = 8)

# 读取seurat_obj  ----                              
seurat_obj <- readRDS("seurat_anno.rds")

gc()


# 小提琴图显示标记基因在不同聚类中的表达----
VlnPlot(seurat_clean, features = c("CD4"), group.by = "cell_type")
RidgePlot(seurat_obj, features = "ESR1", group.by = "phenotype")

# 设定需要查看的 marker genes
genes <- c("LYZ","ITGAE","HOPX",
           "TP63", #基底细胞
           "SCGB1A1", #Clara
           "SFTPD", #AT2
           "FOXJ1", #纤毛细胞
           "AGER", #AT1
           "CSF1R","CD68",
           "IL1B","CD86",
           "MRC1",
           "MARCO","CD163","FABP4","CCR2", #AM
           "CD14","LST1",
           "SIGLEC8","ITGAX",
           "PTPRB","PECAM1",
           "KRT5","KRT14",
           "CD63","KIT",
           "FCGR3B","ITGAM",
           "FCER1A","CST3","XCR1","CLEC9A",
           "TCF4", #pDC
           "GZMB","CD19",
           "MS4A1","CD79A",
           "MZB1","IGKC","JCHAIN",
           # "MKI67","TOP2A","STMN1",
           "GNLY","NKG7","KLRD1",
           "CD3D","CD3E",
           "CD4",
           "CD8A")
genes <- c("TRDC","CTSW","GNLY",
           "LYZ",
           "EPCAM",
           "MS4A1",
           "CPA3",
           "LILRA4","JCHAIN",
           "CD3D","CD8A",
           "CD4","IL7R")
genes <- c("CEPT1")
# 绘制 DotPlot----

p <- DotPlot(
  seurat_obj,
  features = genes,
  assay = "RNA",
  group.by = "condition",
  cols = c("lightgrey", "red"),  # 直接通过 DotPlot 参数设置颜色
  dot.scale = 6                  # 控制最大气泡大小（默认是 6）
) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 8,
      colour = "grey0"
    ),
    axis.text.y = element_text(
      size = 8,
      colour = "grey0"
    ),
    title = element_text(size = 8, colour = "grey0"),
    axis.line = element_line(linewidth = 0.35),
    axis.ticks = element_line(linewidth = 0.35, colour = "grey0"),
    text = element_text(size = 8, colour = "grey0"),  # 全局字体
    axis.title = element_text(size = 8, colour = "grey0")
  )
print(p)
ggsave(filename = "raw_dot_plot.png", plot = p, width = 10, height = 5, dpi = 300)
ggsave(filename = "raw_dot_plot.pdf", plot = p, width = 10, height = 5, dpi = 300)

##分组分亚群绘制----

boxplot_by_celltype <- function(seurat_obj,
                                gene = "MERTK",
                                celltype_col = "cell_type",
                                group_col = "condition",
                                group_levels = c("Control", "Asthma"),
                                color = c("Control" = "lightblue", "Asthma" = "tomato"),
                                output_dir = "output",
                                output_prefix = "MERTK_boxplot",
                                return_label_data = FALSE) {
  
  # 输出路径
  figure_dir <- file.path(output_dir, "annotated_seurat")
  dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 检查基因名是否合法
  if (!(gene %in% rownames(seurat_obj))) {
    stop(paste("Gene", gene, "not found in Seurat object."))
  }
  
  # 提取并整理数据
  data <- FetchData(seurat_obj, vars = c(gene, celltype_col, group_col)) %>%
    mutate(
      !!group_col := factor(.data[[group_col]], levels = group_levels),
      expression = .data[[gene]]
    )
  
  # 计算阳性表达比例
  label_data <- data %>%
    group_by(across(all_of(c(celltype_col, group_col)))) %>%
    summarise(rate = mean(expression > 0) * 100, .groups = "drop")
  
  # 绘图
  p <- ggplot(data, aes(x = .data[[celltype_col]], y = expression, fill = .data[[group_col]])) +
    geom_violin(trim = T, alpha = 0.3, width = 0.8, scale = "width", linewidth = 0.35,
                position = position_dodge(0.8)) +
    # geom_boxplot(outlier.shape = NA, color = "black", size = 0.35, alpha = 0.7,
    #              width = 0.5, position = position_dodge(0.8)) +
    geom_jitter(aes(color = .data[[group_col]]), size = 0.8, alpha = 0.6,
                position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.8)) +
    # stat_summary(fun = mean, geom = "point", color = "blue", size = 2, position = position_dodge(0.8)) +
    # stat_summary(fun = median, geom = "point", color = "red", size = 2, position = position_dodge(0.8)) +
    geom_text(data = label_data, aes(x = .data[[celltype_col]], y = -0.5,
                                     label = paste0(round(rate, 1), "%"),
                                     group = .data[[group_col]], color = .data[[group_col]]),
              position = position_dodge(width = 0.8),
              size = 3, show.legend = FALSE) +
    scale_fill_manual(values = color, guide = "none") +
    scale_color_manual(values = color) +
    labs(x = "Cell Subtype",
         y = bquote(italic(.(gene)) ~ "Relative Expression")) +
    theme_bw() +
    theme(
      legend.position = "top",
      legend.title = element_text(size = 10, face = "bold", color = "black"),
      axis.title = element_text(size = 10, face = "bold", colour = 'black'),
      axis.title.x =  element_blank(),
      axis.text.y = element_text(size = 8, colour = 'black'),
      axis.text.x = element_text(size = 8, angle = 45, hjust = 1, colour = 'black'),
      plot.title = element_text(size = 10, colour = 'black'),
      panel.grid = element_blank(),
      panel.grid.major = element_line(linewidth = 0.35, colour = "grey95")
    )
  
  # 保存图像
  ggsave(file.path(figure_dir, paste0(output_prefix, ".png")), p, width = 6, height = 4, dpi = 300)
  ggsave(file.path(figure_dir, paste0(output_prefix, ".pdf")), p, width = 6, height = 4, dpi = 300)
  
  # 返回图像或阳性比例数据
  if (return_label_data) {
    return(label_data)
  } else {
    return(p)
  }
}


# 绘图----
genelist <- c("S100A9","S100A8")

VlnPlot(MNPs_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "seurat_clusters")



plots <- lapply(genelist, function(g) {
  boxplot_by_celltype(
    seurat_obj = seurat_obj,
    gene = g,
    celltype_col = "seurat_clusters",
    group_col = "condition",
    group_levels = c("Control", "Asthma"),
    color = c("Control" = "lightblue", "Asthma" = "tomato"),
    output_dir = "GCH",
    output_prefix = paste0(g, "_boxplot")
  )
})
plots
# 分细胞类型箱线图----
boxplot_by_celltype_simple <- function(seurat_obj,
                                       gene = "MERTK",
                                       celltype_col = "cell_type",
                                       color = "steelblue",
                                       output_dir = "output",
                                       output_prefix = "MERTK_boxplot",
                                       return_label_data = FALSE) {
  
  # 输出路径
  figure_dir <- file.path(output_dir, "annotated_seurat")
  dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 检查基因是否存在
  if (!(gene %in% rownames(seurat_obj))) {
    stop(paste("Gene", gene, "not found in Seurat object."))
  }
  
  # 提取并整理数据
  data <- FetchData(seurat_obj, vars = c(gene, celltype_col)) %>%
    mutate(expression = .data[[gene]])
  
  # 计算阳性表达比例
  label_data <- data %>%
    group_by(across(all_of(celltype_col))) %>%
    summarise(rate = mean(expression > 0) * 100, .groups = "drop")
  
  # 绘图
  p <- ggplot(data, aes(x = .data[[celltype_col]], y = expression, fill = .data[[celltype_col]])) +
    geom_violin(trim = TRUE, alpha = 0.4, width = 0.8, scale = "width", linewidth = 0.35) +
    geom_jitter(size = 0.8, alpha = 0.6, width = 0.2, color = color) +
    geom_text(data = label_data, aes(x = .data[[celltype_col]], y = -0.5,
                                     label = paste0(round(rate, 1), "%")),
              size = 3, color = "black") +
    scale_fill_manual(values = rep(color, length(unique(data[[celltype_col]])))) +
    labs(x = "Cell Subtype",
         y = bquote(italic(.(gene)) ~ "Relative Expression")) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 10, face = "bold", colour = 'black'),
      axis.title.x =  element_blank(),
      axis.text.y = element_text(size = 8, colour = 'black'),
      axis.text.x = element_text(size = 8, angle = 45, hjust = 1, colour = 'black'),
      plot.title = element_text(size = 10, colour = 'black'),
      panel.grid = element_blank(),
      panel.grid.major = element_line(linewidth = 0.35, colour = "grey95")
    )
  
  # 保存图像
  ggsave(file.path(figure_dir, paste0(output_prefix, ".png")), p, width = 6, height = 4, dpi = 300)
  ggsave(file.path(figure_dir, paste0(output_prefix, ".pdf")), p, width = 6, height = 4, dpi = 300)
  
  # 返回结果
  if (return_label_data) {
    return(label_data)
  } else {
    return(p)
  }
}

FeaturePlot(seurat_obj, features = c("MALT1"))  # 典型浆细胞标志

{
  # AT2标记基因的空间分布
  p1 <- FeaturePlot(seurat_obj, features = c("IL7R"))
  p1
  # p1 <- FeaturePlot(seurat_obj, features = c("SFTA2","SFTPA2", "SFTPD", "NAPSA", "SFTPB", "LAMP3","ABCA3","SFTPC","CLDN18"))
  p1 <- FeaturePlot(seurat_obj, features = c("CD86","EPCAM","PIEZO1","MRC1","MSR1","C1QA"))
  
  p1
  # 保存为 PNG
  ggsave(filename = "FeaturePlot_output.png", plot = p1, width = 8, height = 10, dpi = 300)
  
  # 保存为 PDF
  ggsave(filename = "FeaturePlot_output.pdf", plot = p1, width = 8, height = 10)
}


library(dplyr)
# 按 cluster 和 log fold change 排序，并获取每个 cluster 的前 5 个基因
top5_markers <- markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))

# 查看结果
print(top5_markers)

# 查看前 5 个基因
top5_markers <- head(markers[order(markers$avg_log2FC, decreasing = TRUE), ], 5)
print(top5_markers)

# SingleR,celldex自动注释----
if(T){
  #   # 方法一
  # if (!requireNamespAEC("BiocManager", quietly = TRUE))
  #   install.packages("BiocManager")
  # BiocManager::install("celldex")
  library(SingleR)
  library(celldex)  # 提供参考数据集
  
  # 提取参考数据集
  ref <- celldex::HumanPrimaryCellAtlasData()
  
  # 提取 seurat_obj 数据矩阵
  seurat_obj_matrix <- GetAssayData(AEC_obj, layer = "data")
  
  # 获取聚类信息
  clusters <- AEC_obj$seurat_obj_clusters
  
  # SingleR 注释
  singleR_annotations <- SingleR(test = seurat_obj_matrix, 
                         ref = ref, 
                         labels = ref$label.main,
                         clusters = clusters)
  
  # 将注释结果写入 meta.data
  AEC_obj$SingleR_label <- singleR_annotations$labels[match(AEC_obj$seurat_obj_clusters,
                                                         rownames(singleR_annotations))]
  table(AEC_obj$SingleR_label)
  
  # 方法二
  library(scCATCH)
  # 提取表达矩阵（稀疏矩阵）
  counts_matrix <- GetAssayData(seurat_obj, slot = "counts")
  
  # 提取 cluster 信息
  clusters <- as.character(Idents(seurat_obj))
  
  # 转换为 scCATCH 格式
  data <- createscCATCH(data = counts_matrix, cluster = clusters)
  
  # 使用系统数据库 cellmatch 进行细胞类型注释
  annotations <- findmarkergene(
    data = data,
    species = "Human",
    tissue = "lung",  # 指定组织类型
    marker = "cellmatch"  # 使用系统标记基因数据库
  )
  # 查看注释结果
  print(annotations)
  
}

# 手动定义细胞类型----
# JCHAIN,GZMB, plasma cells
# CD19,MS4A1, B cells
# LYZ MNPs(mononuclear phagocyte)
# CD3D,CD8A, CD8 T cells
# CD4 T cells (CD3D, CD4, and IL7R)
# AECs (EPCAM,KRT78)
# natural killer cells (GNLY)
# mast cells (CPA3)
if(T){
  
  
  # 初始化所有 cluster 的注释为原始编号
  cluster_annotations <- as.character(Idents(seurat_obj))
  
  # 手动定义部分 cluster 的细胞类型
  {cluster_annotations[Idents(seurat_obj) %in% c(3,4,5,7,15)] <- "AECs"  
  cluster_annotations[Idents(seurat_obj) %in% c(6,9,10,11,12,14)] <- "MNPs"   
  cluster_annotations[Idents(seurat_obj) %in% c(0,1)] <- "CD8 T cells"
  cluster_annotations[Idents(seurat_obj) %in% 2] <- "CD4 T cells"
  cluster_annotations[Idents(seurat_obj) %in% 8] <- "B cells"
  # cluster_annotations[Idents(seurat_obj) %in% 11] <- "plasma cells"
  cluster_annotations[Idents(seurat_obj) %in% 13] <- "Mast cells"
  cluster_annotations[Idents(seurat_obj) %in% 16] <- "NK cells"
  }
  # 其他 cluster 暂不定义，保留原始编号
  
  # 分配细胞类型
  seurat_obj$cell_type <- cluster_annotations
  
  library(RColorBrewer)
  
  # 生成颜色
  cell_types <- unique(seurat_obj$cell_set_clean)
  colors <- colorRampPalette(brewer.pal(10, "Set3"))(length(cell_types))
  
  p <- DimPlot(seurat_obj, reduction = "umap", group.by = "cell_set_clean", order = T, label = T, label.size = 3) +
    labs(title = "Global clustering") +
    scale_color_manual(values = setNames(colors, cell_types)) +
  theme_minimal() +  #  theme_classic()
    theme(
      text = element_text(size = 8, face = "bold", hjust = 0.5),
      plot.title = element_text(size = 8, face = "bold", hjust = 0.5),  # 标题居中 & 加粗
      axis.title = element_blank(),  # 隐藏坐标轴标题
      axis.text = element_blank(),
      panel.grid = element_blank(),  # 去除网格线
      legend.position = "right"      # 图例放右侧（可改 "bottom", "top"）
    )
  # # 用 geom_encircle 添加包围圈
  # p + geom_encircle(aes(x = UMAP_1, y = UMAP_2, group = cell_type), 
  #                   data = umap_df, 
  #                   expand = 0.05,  # 控制边界宽度
  #                   color = "black", 
  #                   linetype = "dashed", 
  #                   size = 1)
  p
  ggsave(filename = "UMAP_All.png", plot = p, width = 6, height = 5, dpi = 300)
  ggsave(filename = "UMAP_All.pdf", plot = p, width = 6, height = 5, dpi = 300)
}

getwd()
# 分组查看基因表达----
# 提取AA组的数据
seurat_Ag <- subset(seurat, subset = sample == "Ag")

seurat_Pre <- subset(seurat, subset = sample == "Pre")
seurat_Dil <- subset(seurat, subset = sample == "Dil")
# 提取ANA组的数据
seurat_AA <- subset(seurat, subset = phenotype == "AA")
seurat_ANA <- subset(seurat, subset = phenotype == "ANA")
# 设置分组身份，差异分析将默认此分组
Idents(seurat) <- "phenotype"
# 可视化AA组的Gene表达
gene <- "MALT1"
p <- FeaturePlot(seurat_obj, 
            features = gene, 
            # split.by = "condition",
            reduction = "umap", 
            cols = c("lightgrey", "red"),  # 低表达灰色，高表达红色
            pt.size = 0.1,
            order = TRUE)   # 按表达量排序，高表达细胞覆盖在顶部
  # ggtitle(gene)
p 
ggsave(
  filename = "CEPT1_AEC_FeaturePlot.png",  
  plot = p,               
  width = 5,                       
  height = 5,                        
  dpi = 300                          
)
ggsave(
  filename = "CEPT1_AEC_FeaturePlot.pdf",  
  plot = p,               
  width = 5,                       
  height = 5,                        
  dpi = 300                          
)

# 可视化ANA组的Gene表达
FeaturePlot(seurat_ANA, 
            features = "IL1B", 
            reduction = "umap", 
            cols = c("lightgrey", "blue"),  # 颜色可自定义区分组别
            pt.size = 1.5, 
            order = TRUE) +
  ggtitle("ANA组 - GeneX表达")

# 不同组间同一基因的表达
theme_set(theme_bw(base_size = 8)) 
p <- FeaturePlot(seurat, 
            features = "TSSK4", 
            split.by = "phenotype",  # 按group分组分面
            cols = c("grey90", "red"),
            pt.size = 0.5,
            order = TRUE,
            ncol = 2  # 分2列排列
            ) 
p 
ggsave(
  filename = "CSF2_FeaturePlot.png",  
  plot = p,               
  width = 6,                       
  height = 4,                        
  dpi = 300                          
)
ggsave(
  filename = "CSF2_FeaturePlot.pdf",  
  plot = p,               
  width = 6,                       
  height = 4,                        
  dpi = 300                          
)

# 在metadata中标记AA组
seurat@meta.data$is_AA <- ifelse(seurat$phenotype == "AA", 1, 0)
# 可视化AA组位置及GeneX表达叠加效果
FeaturePlot(seurat, features = c("is_AA", "IL1B"), 
            cols = c("grey", "#FF0000", "blue"),
            blend = T)


# 按聚类查看差异基因----
# 查看当前聚类标识
head(Idents(seurat)) 
# 若未设置，手动指定聚类结果（例如 "seurat_clusters"）
Idents(seurat) <- seurat$seurat_clusters
# 查找 Cluster 0 相对于 Cluster 1 上调的基因（默认方法：Wilcoxon 秩和检验）
diff_genes <- FindMarkers(
  seurat,
  ident.1 = 8,          # 目标亚类（例如 Cluster 0）
  ident.2 = 7,          # 对照亚类（例如 Cluster 1）
  logfc.threshold = 0.25,  # 过滤最小 log2 倍变化（根据数据调整）
  min.pct = 0.1,        # 基因在至少 10% 的细胞中表达
  test.use = "wilcox"   # 使用 Wilcoxon 秩和检验（默认）
)

# 按调整后 p 值排序（最显著在前）
diff_genes <- diff_genes[order(diff_genes$p_val_adj), ]

# 提取 log2 倍变化（avg_log2FC）绝对值最大的前 10 个基因
top_genes <- diff_genes[order(abs(diff_genes$avg_log2FC), decreasing = TRUE), ]
head(top_genes, 10)

# 提取前 20 个差异基因
top20 <- rownames(diff_genes)[1:20]

# 绘制热图
DoHeatmap(
  seurat,
  features = top20,
  group.by = "seurat_clusters",
  slot = "scale.data"  # 使用标准化后的数据
)

FeaturePlot(
  seurat_obj,
  features = "ESR1",  # 前 4 个基因
  # split.by = "condition",
  cols = c("grey90", "red"),
  order = TRUE
)


# 先筛选出感兴趣的亚群（例如 cluster0 和 cluster1）
seurat_sub <- subset(seurat, seurat_clusters %in% c(2,3,5,12))  # 替换成你需要的cluster编号

# 然后对子集数据画图，split.by参数会基于筛选后的亚群分割
FeaturePlot(
  AEC_obj,
  features = c("CD3E","CD8A"), 
  # split.by = "seurat_clusters",  # 此时只会显示c(0,1)两个亚群
  # idents = c(4, 7),
  cols = c("grey90", "red"),
  order = TRUE
)

VlnPlot(
  AEC_obj,
  features = c("CD3E","CD8A"),  # 例如第一个基因
  pt.size = 0.5,
  # idents = c(4, 7)      # 仅显示目标亚类
)



if(T){
  # 绘制细胞类型的条形图或饼图，展示各类型的比例----
  # 动态生成所需数量的颜色
  # 生成颜色
  cell_types <- unique(seurat$cell_set_clean)
  colors <- colorRampPalette(brewer.pal(6, "Set2"))(length(cell_types))
  
  # 应用到 DimPlot
  p2 <- DimPlot(seurat, reduction = "umap", group.by = "cell_set_clean", label = TRUE, label.size = 4) +
    labs(title = "UMAP with Annotated Cell Types") +
    scale_color_manual(values = setNames(colors, cell_types)) +
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      # title = element_blank(),
      axis.line = element_line(linewidth = 0.35),
      axis.ticks = element_line(linewidth = 0.35, colour = "grey0"),
      text = element_text(size = 8, colour = "grey0"),  # 设置全局字体大小
      axis.text.x = element_text(size = 8, colour = "grey0"),  # 明确设置坐标轴刻度标签字体大小
      axis.text.y = element_text(size = 8, colour = "grey0"),  # 明确设置坐标轴刻度标签字体大小
      axis.title = element_text(size = 8, colour = "grey0"),  # 明确设置坐标轴标题字体大小
      legend.position = "right") 
  p2
  
  # 保存
  ggsave(paste0("细胞类型", ".tiff"), plot = p2, height = 6, width = 6,dpi = 300)
  ggsave(paste0("细胞类型", ".pdf"), plot = p2, height = 6, width = 6)
  ##饼图----
  df <- as.data.frame(table(seurat$cell_type))
  p3 <- ggplot(df, aes(x = "", y = Freq, fill = Var1)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y") +
    scale_fill_manual(values = setNames(colors, unique(df$Var1))) +
    labs(title = "Cell Type Proportion") +
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      title = element_text(size = 8, colour = "grey0"),
      # axis.line = element_line(linewidth = 0.35),
      axis.ticks = element_line(linewidth = 0.35, colour = "grey0"),
      text = element_text(size = 8, colour = "grey0"),  # 设置全局字体大小
      axis.text.x = element_text(size = 8, colour = "grey0"),  # 明确设置坐标轴刻度标签字体大小
      axis.text.y = element_text(size = 8, colour = "grey0"),  # 明确设置坐标轴刻度标签字体大小
      axis.title = element_blank(),  # 明确设置坐标轴标题字体大小
      legend.position = "right") 
  p3
  ggsave(paste0("细胞类型比例", ".tiff"), plot = p3, height = 8, width = 8,dpi = 300)
  ggsave(paste0("细胞类型比例", ".pdf"), plot = p3, height = 8, width = 8)
  # 条形图----
  if(T){
    # cell proportion
    df <- AEC_obj@meta.data
    df_prop <- df %>%
      dplyr::group_by(cell_types, condition, Channel) %>%
      dplyr::summarise(Count = dplyr::n(), .groups = "drop") %>% 
      as.data.frame() %>% 
      setNames(c("Lineage", "condition", "sample","Count"))
    
    df_prop <- df_prop %>%
      mutate(
        condition = factor(condition, levels = c("Control", "Asthma")),
        sample = factor(sample, levels = unique(sample[order(condition)]))  # 按condition排序
      )
    df_prop <- arrange(df_prop, condition)
    
    p1 <-  df_prop %>%
      ggplot(aes(y = forcats::fct_rev(forcats::fct_infreq(Lineage)), 
                 x = Count,
                 fill = Lineage)) +
      geom_col() +
      # geom_bar(stat = 'Count') +
      labs(x = 'Cell count', y = NULL) +
      scale_fill_manual(
        name = "Cell lineage",
        values = colors,
        labels = names(colors)
      ) +
      theme_bw(base_size = 10) +
      theme(
        axis.text = element_text(size = 8, color = 'black'),
        legend.position = "none"  # Remove legend as it will be in the right panel
      )

    # Create frequency plot (right panel)
    p2 <- ggplot(df_prop, 
                 aes(x = Count,
                     y = condition,
                     fill = Lineage)) +
      geom_bar(position = "fill", stat = "identity") +
      # coord_flip() +
      labs(y = NULL, x = 'Cell Lineage Proportion') +
      # scale_x_discrete(labels = sample_labels) +
      scale_fill_manual(
        name = "Cell lineage",
        values = colors
      ) +
      theme_bw(base_size = 10) +
      theme(
        axis.text = element_text(size = 8, color = 'black'),
        legend.key.size = unit(0.5, "cm")  # Adjust legend key size
      ) 
    library(patchwork)
    # Combine plots
    coplot <- (p1 + p2) + 
      plot_layout(guides = "collect") & 
      plot_annotation(
        title = "Proportion of cell classification",
        theme = theme(
          plot.title = element_text(
            hjust = 0.5, 
            size = 10, 
            face = 'bold'
          )
        )
      )
  }
  coplot
  
  
  # Save plot
  ggsave(
    filename = "AEC分组细胞比例图2.pdf",
    plot = coplot,
    width = 8,  # Slightly wider to accommodate legend
    height = 5   # Adjusted aspect ratio
  )
  ggsave(
    filename = "AEC分组细胞比例图2.png",
    plot = coplot,
    width = 8,  # Slightly wider to accommodate legend
    height = 5   # Adjusted aspect ratio
  )

  
# 备份----
saveRDS(seurat, file = "seurat_Dim.rd")

seurat <- readRDS("seurat_Dim.rds")
str(seurat)
# 1. 提取 UMAP 坐标和关键元数据----
plot_data <- FetchData(MNPs, vars = c("UMAP_1", "UMAP_2", "cell_type", "condition", "phenotype"))
# 安装 ggpointdensity（如果未安装）
if (!require(ggpointdensity)) install.packages("densmap")

# 绘制点密度图（显示细胞分布）
p_density <- ggplot(MNPs, aes(x = UMAP_1, y = UMAP_2)) +
  geom_pointdensity(size = 0.5, alpha = 0.5) +  # 点密度
  scale_color_viridis_c(name = "Density") +     # 密度颜色
  facet_grid(
    rows = vars(disease),
    cols = vars(group)
  ) +
  theme_bw()

# 计算各 MNP 亚群比例
prop_data <- umap_data %>%
  group_by(disease, group, MNP_subset) %>%
  summarise(n = n()) %>%
  mutate(prop = n / sum(n))  # 计算比例

# 绘制比例条形图
p_prop <- ggplot(prop_data, aes(x = MNP_subset, y = prop, fill = MNP_subset)) +
  geom_bar(stat = "identity") +
  facet_grid(
    rows = vars(disease),
    cols = vars(group)
  ) +
  labs(y = "Proportion", x = "MNP Subset") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 组合图形
combined_plot <- p_density / p_prop + plot_layout(heights = c(2, 1))
print(combined_plot)

# 提取上皮细胞亚群再细分----
library(Harmony)
if(T){
  AEC <- subset(seurat_obj, subset = cell_set_clean == "Epithelial cells")
  # 查看所有的meta.data列名
  colnames(AEC@meta.data)
  # 数据标准化处理流程+细胞注释----
  
  # marker_list <- list(
  #   AT1 = c("AGER"),
  #   AT2 = c("SFTPD"),
  #   Ciliated_Cells = c("FOXJ1","PIFO"),
  #   Goblet_cells = c("MUC5AC"),
  #   Basal_cells = c("KRT5", "S100A2"),
  #   Clara_cells = c("SCGB1A1", "GNLY")
  # )
  
  subType <- function(x, regress_cell_cycle = TRUE, marker_list = NULL) {
    # Step 1: Normalize
    x <- NormalizeData(x, normalization.method = "LogNormalize")
    
    # Step 2: Find variable features
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
    
    # Step 3: Calculate cell cycle scores,optional
    if (regress_cell_cycle) {
      s_genes <- cc.genes$s.genes
      g2m_genes <- cc.genes$g2m.genes
      x <- CellCycleScoring(x, s.features = s_genes, g2m.features = g2m_genes)
      
      # Step 4: Scale and regress out nCount_RNA and cell cycle
      x <- ScaleData(x, vars.to.regress = c("nCount_RNA", "S.Score", "G2M.Score"))
    } else {
      x <- ScaleData(x, vars.to.regress = "nCount_RNA")
    }
    
    # Step 5: PCA
    x <- RunPCA(x)
    
    # Step 6: Harmony batch correction
    # x <- RunHarmony(x, group.by.vars = "Dataset")
    
    # Step 7: UMAP
    x <- RunUMAP(x,
                 # reduction = "harmony", 
                 dims = 1:20, seed.use = 12345)
    # ElbowPlot(MNPs)
    # Step 8: Clustering
    x <- FindNeighbors(x, 
                       # reduction = "harmony",
                       dims = 1:20, verbose = FALSE)
    x <- FindClusters(x, resolution = 0.3, verbose = FALSE, random.seed = 20220727)
    return(x)
  }
  
  AEC_obj <- subType(AEC)
  table(AEC_obj$seurat_clusters)
  DimPlot(AEC_obj, group.by = "seurat_clusters", label = TRUE)

  if (!is.null(marker_list)) {
    for (cell_type in names(marker_list)) {
      FeaturePlot(AEC_obj, features = marker_list[[cell_type]], label = TRUE) +
        ggtitle(paste("Marker genes for", cell_type))
    }
  }
  
  AEC_obj_markers <- FindAllMarkers(AEC_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  write.csv(AEC_obj_markers,file = "AEC_obj_markers.csv", row.names = T)
  head(AEC_obj_markers)
  head(AEC_obj)
  
  # 生成颜色
  table(AEC_obj$Channel)
  seurat_clusters <- unique(AEC_obj$seurat_clusters)

  colors <- colorRampPalette(brewer.pal(8, "Set3"))(length(seurat_clusters))
  # 自定义配色
 p <- DimPlot(AEC_obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE, label.size = 4) +
    labs(title = "UMAP OF MNPs") +
    scale_color_manual(values = setNames(colors, seurat_clusters))
 p
 ggsave(filename = "AEC.png",plot = p, width = 5.5, height = 5) 
 ggsave(filename = "AEC.pdf",plot = p, width = 5.5, height = 5) 
  # # 差异基因分析
  at2_markers <- FindMarkers(MNPs, ident.1 = "AT2 cell", ident.2 = NULL, group.by = "cell_type")
  # 查看结果
  head(at2_markers)
  table(MNPs$phenotype)

  
  if(T){
    


    # # 手动定义部分 cluster 的细胞类型,上皮
    genes <- c("TP63", #基底细胞
               "TYMS","CENPF","MKI67","KRT5",
               "CEPT1",
               "LYZ","LTF","MUC5B",
               "SCGB1A1","SCGB3A1", #Clara
               "SFTPD", #AT2
               "FOXJ1","PIFO", #Ciliated Cells
               "AGER","AGR3","C2orf40",
               
               "MUC21","MUC5AC" #Goblet cells
    )


    # # 手动定义部分 cluster 的细胞类型,巨噬
    if(T){
      genes <- c("PIEZO1","ITGAM",
                 
                 "CEACAM8", "FCGR3A", #CD66B(CEACAM8),NEU
                 # "CD19","CD20","CD38","CD138","CD3E","CD4","CD8A",
                 "IL1B","CD86","CD163",
                 "IRF8", "BATF3",
                 "TCF4","GZMB","JCHAIN","IL3RA", #pDC
                 "CD1C","CD1E","ITGAX","CD1B","FCER1A", "CST3", #cDC
                 
                 "SERPINF1",
                 "MRC1",
                 "CSF1R","CD68",
                 "CD14","FCN1",
                 "MARCO","MSR1","C1QA","C1QB",
                 "CCL18","HOPX","PPARG","CCR2", #AM
                 "LYZ","FABP4","APOE","VEGFA","AREG","MERTK","CCL2","CCL3"
      )
    }


    # 绘制 DotPlot----
    p <- DotPlot(AEC_obj, 
                 features = "CEPT1", 
                 group.by = "seurat_clusters",
                 scale = T,
                 scale.by = "size",
                 col.min = -2,
                 col.max = 3,
                 dot.min = 0,
                 cols = c("grey", "red"),
                 dot.scale = 5) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        title = element_text(size = 8, colour = "grey0"),
        axis.line = element_line(linewidth = 0.35),
        axis.ticks = element_line(linewidth = 0.35, colour = "grey0"),
        text = element_text(size = 8, colour = "grey0"),  # 设置全局字体大小
        axis.text = element_text(size = 8, colour = "grey0"),  # 明确设置坐标轴刻度标签字体大小
        axis.title = element_text(size = 8, colour = "grey0"))
    p
    ggsave(filename = "AEC_dot_plot.png", plot = p, width = 10, height = 5, dpi = 300)
    ggsave(filename = "AEC_dot_plot.pdf", plot = p, width = 10, height = 5, dpi = 300)
    
    # 初始化所有 cluster 的注释为原始编号
    Idents(AEC_obj) <- "seurat_clusters"
    cluster_annotations <- as.character(Idents(AEC_obj))
    
    {cluster_annotations[Idents(AEC_obj) %in% c(1)] <- "Clara Cells"
      cluster_annotations[Idents(AEC_obj) %in% c(2)] <- "Goblet cells"
      cluster_annotations[Idents(AEC_obj) %in% c(3)] <- "Serous Cells"
      cluster_annotations[Idents(AEC_obj) %in% c(4)] <- "Multiciliated cell"
      cluster_annotations[Idents(AEC_obj) %in% c(5)] <- "Mucous-ciliated"
      cluster_annotations[Idents(AEC_obj) %in% c(0,7)] <- "Basal cells"
      cluster_annotations[Idents(AEC_obj) %in% 6] <- "Cycling basal cells"
      cluster_annotations[Idents(AEC_obj) %in% 8] <- "Squamous metaplasia"
      cluster_annotations[Idents(AEC_obj) %in% 9] <- "Ciliated Cells"
      cluster_annotations[Idents(AEC_obj) %in% 10] <- "Ionocytes"
      cluster_annotations[Idents(AEC_obj) %in% 11] <- "Precursor ciliated cells"
    }
    # 分配细胞类型
    AEC_obj$cell_types <- cluster_annotations
    table(AEC_obj$cell_types)
    
    library(RColorBrewer)
    
    # 生成颜色
    cell_types <- unique(AEC_obj$cell_types)
    colors <- colorRampPalette(brewer.pal(11, "Set3"))(length(cell_types))
    # 自定义配色
    p <- DimPlot(AEC_obj, reduction = "umap", group.by = "cell_types", 
                 pt.size = 0.01, alpha = 1,
                 order = T, label = T, label.size = 3) +
      labs(title = "MNPs clustering") +
      scale_color_manual(values = setNames(colors, cell_types)) +
      theme_minimal() +  #  theme_classic()
      theme(
        text = element_text(size = 8, face = "bold", hjust = 0.5),
        plot.title = element_text(size = 8, face = "bold", hjust = 0.5),  # 标题居中 & 加粗
        axis.title = element_blank(),  # 隐藏坐标轴标题
        axis.text = element_blank(),
        panel.grid = element_blank(),  # 去除网格线
        legend.position = "none"      # 图例放右侧（可改 "bottom", "top"）
      )
    p
    ggsave(filename = "UMAP_AEC_obj.png", plot = p, width = 5, height = 5, dpi = 300)
    ggsave(filename = "UMAP_AEC_obj.pdf", plot = p, width = 5, height = 5, dpi = 300)
  }
  
  VlnPlot(AEC_obj, features = "CEPT1", group.by = "cell_types")
  AEC_ciliated_obj <- subset(AEC_obj, cell_types %in% c("Multiciliated cell"))
  # 绘制 DotPlot----
  p <- DotPlot(AEC_ciliated_obj, 
               features = "CEPT1", 
               # split.by = "condition",
               group.by = "condition",
               # scale = T,
               # scale.by = "size",
               col.min = -2,
               col.max = 3,
               dot.min = 0,
               cols = c("grey", "red"),
               dot.scale = 5) +
    
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      title = element_text(size = 8, colour = "grey0"),
      axis.line = element_line(linewidth = 0.35),
      axis.ticks = element_line(linewidth = 0.35, colour = "grey0"),
      text = element_text(size = 8, colour = "grey0"),  # 设置全局字体大小
      axis.text = element_text(size = 8, colour = "grey0"),  # 明确设置坐标轴刻度标签字体大小
      axis.title = element_text(size = 8, colour = "grey0"))
  p
  ggsave(filename = "raw_dot_plot.tiff", plot = p, width = 10, height = 5, dpi = 300)
  ggsave(filename = "raw_dot_plot.pdf", plot = p, width = 10, height = 5, dpi = 300)
  
  
  saveRDS(AEC_obj,"AEC_obj.rds")
  
  ##分组分亚群绘制----
  #提取数据（必须包含cell_type和condition）
  data <- FetchData(AEC_obj, vars = c("CEPT1", "cell_types", "condition"))
  data <- data %>% mutate(condition = factor(condition, levels = c("Control", "Asthma")))
  data <- arrange(data, condition)
  # 按细胞亚类和组别计算表达比例
  label_data <- data %>%
    group_by(cell_types, condition) %>%
    summarise(rate = mean(CEPT1 > 0) * 100, .groups = "drop")
  
  p <- ggplot(data, aes(x = cell_types, y = CEPT1, fill = condition)) +
    geom_boxplot(
      outlier.shape = NA,
      color = "black",
      size = 0.35,
      alpha = 0.7,
      width = 0.5,
      position = position_dodge(width = 0.9)
    ) +
    geom_violin(trim = T, alpha = 0.3, width = 0.9, scale = "width", linewidth = 0.35) + # 增加小提琴图，透明度调低
    geom_jitter(
      aes(color = condition),
      size = 1,
      alpha = 0.6,
      position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.9)
    ) +
    
    stat_summary(fun = mean, geom = "point", color = "blue", size = 2, position = position_dodge(0.9)) +
    stat_summary(fun = median, geom = "point", color = "red", size = 2, position = position_dodge(0.9)) +
    # stat_compare_means(
    #   aes(group = condition),
    #   label = "p.format",
    #   size = 3,
    #   method = "kruskal.test"
    # ) +
    geom_text(
      data = label_data,
      aes(x = cell_types, y = -0.5, label = paste0(round(rate, 1), "%"), group = condition, color = condition),
      position = position_dodge(width = 0.9),
      size = 3,
      show.legend = FALSE
    ) +
    scale_fill_manual(values = c("Control" = "lightblue", "Asthma" = "tomato"),guide = "none") +
    scale_color_manual(values = c("Control" = "lightblue", "Asthma" = "tomato")) +
    labs(
      title = "",
      x = "Cell Subtype",
      y = bquote(italic(.("CEPT1")) ~ " Relative Expression")
    ) +
    theme_bw() +
    theme(
      legend.position = "top",
      axis.text.y  = element_text(size = 8,  colour = 'black'),
      axis.text.x = element_text(size = 8, angle = 45,hjust = 1, colour = 'black'),
      panel.grid = element_blank(),
      panel.grid.major = element_line(linewidth = 0.35, colour = "grey95")
    )
  
  p
  
  #分类不分组----
  data <- FetchData(AEC_obj, vars = c("CEPT1", "cell_types", "condition"))
  data <- data %>% mutate(condition = factor(condition, levels = c("Control", "Asthma")))
  data <- arrange(data, condition)
  # 按细胞亚类和组别计算表达比例
  label_data <- data %>%
    group_by(cell_types) %>%
    summarise(rate = mean(CEPT1 > 0) * 100, .groups = "drop")
  
  p <- ggplot(data, aes(x = cell_types, y = CEPT1, fill = cell_types)) +
    geom_boxplot(
      outlier.shape = NA,
      color = "black",
      size = 0.35,
      alpha = 0.7,
      width = 0.5,
      position = position_dodge(width = 0.9)
    ) +
    geom_violin(trim = T, alpha = 0.3, width = 0.9, scale = "width", linewidth = 0.35) + # 增加小提琴图，透明度调低
    geom_jitter(
      aes(color = cell_types),
      color = "black",
      size = 1,
      alpha = 0.5,
      position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.9)
    ) +
    
    stat_summary(fun = mean, geom = "point", color = "blue", size = 2, position = position_dodge(0.9)) +
    stat_summary(fun = median, geom = "point", color = "red", size = 2, position = position_dodge(0.9)) +
    # stat_compare_means(
    #   aes(group = condition),
    #   label = "p.format",
    #   size = 3,
    #   method = "kruskal.test"
    # ) +
    geom_text(
      data = label_data,
      aes(x = cell_types, y = -0.5, label = paste0(round(rate, 1), "%"), group = cell_types, color = cell_types),
      color = "black",
      position = position_dodge(width = 0.9),
      size = 3,
      show.legend = FALSE
    ) +
    scale_fill_manual(values = colors, guide = "none") +
    # scale_color_manual(values = colors, guide = "none") +
    labs(
      title = "",
      x = "Cell Subtype",
      y = bquote(italic(.("CEPT1")) ~ " Relative Expression")
    ) +
    theme_bw() +
    theme(
      legend.position = "top",
      text = element_text(size = 8,  colour = 'black'),
      axis.text.y  = element_text(size = 8,  colour = 'black'),
      axis.text.x = element_text(size = 8, angle = 45,hjust = 1, colour = 'black'),
      panel.grid = element_blank(),
      panel.grid.major = element_line(linewidth = 0.35, colour = "grey95")
    )
  
  p
  
  ggsave(filename = "CEPT1上皮不同细胞只分类占比.png",  width = 11, height = 5, dpi = 300)
  ggsave(filename = "CEPT1上皮不同细胞只分类占比.pdf",  width = 11, height = 5, dpi = 300)
  
  # Load R packages
  library(glue)
  library(dplyr)
  library(NMF)
  library(reshape2)
  library(ComplexHeatmap)
  # Running NMF algorithms, 为了快速复现，这里设置将原来 nrun = 100 设置为 nrun = 1
  ranks <- 2:10
  estim.coad <- nmf(scale_ratio, ranks, nrun = 1, method = "lee")
  
  
  expression_data <- as.data.frame(FetchData(AEC, vars = c("HTRA2", "condition")))
  
  # 检查数据
  head(expression_data)
  
  # 小提琴图
  expression_data <- expression_data %>% mutate(group = factor(condition, levels = c("Healthy", "Disease")))
  expression_data <- arrange(expression_data, group)
  # 移除表达为0的样本
  # expression_data <- expression_data[expression_data$ABLIM1 > 0, ]
  
  # 绘制箱线图
  p <- ggplot(expression_data, aes(x = group, y = HTRA2, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.3, scale = "width", linewidth = 0.35) + # 增加小提琴图，透明度调低
    geom_boxplot(aes(fill = group),
                 alpha = 0.5,
                 linewidth = 0.5,
                 notch = FALSE,
                 width = 0.5,
                 outlier.size = 0.05,
                 outlier.color = "black"
    ) +
    geom_jitter(mapping = aes(colour = group),
                shape = 16,
                position = position_jitter(0.2),
                size = 0.3,
                alpha = 0.8
    ) +
    geom_signif(
      comparisons = list(c("Healthy", "Disease")), 
      map_signif_level = T, 
      test = "t.test", # 或 "wilcox.test" 视具体情况
      step_increase = 1,
      tip_length = 0,
      textsize = 3.5,
      size = 0.35
    ) +
    stat_summary(fun = mean, geom = "point", color = "blue", size = 2) +
    stat_summary(fun = median, geom = "point", color = "red", size = 2) +
    labs(
      title = "HTRA2 Expression in Healthy vs Disease", 
      x = "Condition", 
      y = bquote(italic(.("HTRA2")) ~ "mRNA (fold change)")
    ) +
    theme(
      # panel.grid.major = element_line(color = "grey90", size = 0.35),
      # panel.border = element_rect(color = "black", fill = NA, size = 0.35),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line(linewidth = 0.35),
      axis.ticks = element_line(linewidth = 0.35, colour = "grey0"),
      text = element_text(size = 8, colour = "grey0"),
      axis.text.x = element_text(size = 8,  colour = "grey0"),
      axis.text.y = element_text(size = 8, colour = "grey0"),
      axis.title.y = element_text(size = 8, colour = "grey0"),
      axis.title.x = element_blank(),
      plot.title = element_text(size = 8, fAEC = "bold"),
      legend.position = "none"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
    # scale_y_log10(expand = expansion(mult = c(0, 0.1))) +
    scale_fill_manual(values = c("Healthy" = "#87cEEB", "Disease" = "#FF5722")) +
    scale_colour_manual(values = c("Healthy" = "#87cEEB", "Disease" = "#FF5722"))
  # 显示图
  print(p)
  # 保存图片,多个柱子，宽度加0.5
  ggsave(paste0("total_HTRA2", ".tiff"), plot = p, height = 3, width = 3,dpi = 300)
  ggsave(paste0("total_HTRA2", ".pdf"), plot = p, height = 3, width = 3,dpi = 300)
  
  result <- t.test(HTRA2~group, expression_data, var.equal = T)#方差齐
  print(result)
  p1 <- FeaturePlot(AEC, features = c("HTRA2"))
  p1
}

# 提取AM亚群再细分----
if(T){
  AM <- subset(AEC, subset = cell_type == "AM")
  # 查看提取的结果
  print(AM)
  # 标准化数据
  AM <- NormalizeData(AM)
  
  # 识别高变基因
  AM <- FindVariableFeatures(AM)
  
  # 缩放数据
  AM <- ScaleData(AM)
  
  # PCA 降维
  AM <- RunPCA(AM)
  ElbowPlot(AM)
  # 最近邻图和聚类
  AM <- FindNeighbors(AM, dims = 1:10)
  AM <- FindClusters(AM, resolution = 0.5)
  
  # UMAP 降维
  AM <- RunUMAP(AM, dims = 1:10)
  
  # 可视化新的聚类结果
  DimPlot(AM, reduction = "umap", label = TRUE)
  
  AM_sub_markers <- FindAllMarkers(AM, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  head(AM_sub_markers)
  head(AM)
  
  # 生成颜色
  cell_type <- unique(AM$seurat_clusters)
  colors <- colorRampPalette(brewer.pal(12, "Set3"))(length(cell_type))
  # 自定义配色
  DimPlot(AM, reduction = "umap", group.by = "seurat_clusters", label = TRUE, label.size = 4) +
    labs(title = "UMAP OF AT2") +
    scale_color_manual(values = setNames(colors, cell_type))
  
  # # 差异基因分析
  AM_sub_markers <- FindMarkers(seurat, ident.1 = "AT2 cell", ident.2 = NULL, group.by = "cell_type")
  # 查看结果
  head(AM_sub_markers)
  table(AM$condition)
  
  # 添加分组信息
  AM$condition <- ifelse(grepl("Donor", AM$sample), "Healthy", "Disease")
  table(AM$condition)
  # 设置分组
  Idents(AM) <- "condition"
  
  expression_data <- as.data.frame(FetchData(AM, vars = c("PIEZO1", "condition")))
  
  # 检查数据
  head(expression_data)
  
  # 小提琴图
  expression_data <- expression_data %>% mutate(group = factor(condition, levels = c("Healthy", "Disease")))
  expression_data <- arrange(expression_data, group)
  # 移除表达为0的样本
  # expression_data <- expression_data[expression_data$ABLIM1 > 0, ]
  
  # 绘制箱线图
  p <- ggplot(expression_data, aes(x = group, y = PIEZO1, fill = group)) +
    geom_violin(trim = FALSE, alpha = 0.3, scale = "width", linewidth = 0.35) + # 增加小提琴图，透明度调低
    geom_boxplot(aes(fill = group),
                 alpha = 0.5,
                 linewidth = 0.5,
                 notch = FALSE,
                 width = 0.5,
                 outlier.size = 0.05,
                 outlier.color = "black"
    ) +
    geom_jitter(mapping = aes(colour = group),
                shape = 16,
                position = position_jitter(0.2),
                size = 0.3,
                alpha = 0.8
    ) +
    geom_signif(
      comparisons = list(c("Healthy", "Disease")), 
      map_signif_level = T, 
      test = "t.test", # 或 "wilcox.test" 视具体情况
      step_increase = 1,
      tip_length = 0,
      textsize = 3.5,
      size = 0.35
    ) +
    stat_summary(fun = mean, geom = "point", color = "blue", size = 2) +
    stat_summary(fun = median, geom = "point", color = "red", size = 2) +
    labs(
      title = "PIEZO1 Expression in Healthy vs Disease", 
      x = "Condition", 
      y = bquote(italic(.("PIEZO1")) ~ "mRNA (fold change)")
    ) +
    theme(
      # panel.grid.major = element_line(color = "grey90", size = 0.35),
      # panel.border = element_rect(color = "black", fill = NA, size = 0.35),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line(linewidth = 0.35),
      axis.ticks = element_line(linewidth = 0.35, colour = "grey0"),
      text = element_text(size = 8, colour = "grey0"),
      axis.text.x = element_text(size = 8,  colour = "grey0"),
      axis.text.y = element_text(size = 8, colour = "grey0"),
      axis.title.y = element_text(size = 8, colour = "grey0"),
      axis.title.x = element_blank(),
      plot.title = element_text(size = 8, fAEC = "bold"),
      legend.position = "none"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
    # scale_y_log10(expand = expansion(mult = c(0, 0.1))) +
    scale_fill_manual(values = c("Healthy" = "#87cEEB", "Disease" = "#FF5722")) +
    scale_colour_manual(values = c("Healthy" = "#87cEEB", "Disease" = "#FF5722"))
  # 显示图
  print(p)
  # 保存图片,多个柱子，宽度加0.5
  ggsave(paste0("total_HTRA2", ".tiff"), plot = p, height = 3, width = 3,dpi = 300)
  ggsave(paste0("total_HTRA2", ".pdf"), plot = p, height = 3, width = 3,dpi = 300)
  
  result <- t.test(HTRA2~group, expression_data, var.equal = T)#方差齐
  print(result)
  p1 <- FeaturePlot(AM, features = c("HTRA2"))
  p1
}

# 提取AT细胞亚群再细分----

if(T){
  AECs <- readRDS("AEC.rds")
  # 提取 AT cell
  AECs <- subset(seurat_obj, subset = cell_type == "Airway epithelial cells")
  # 查看提取的结果
  print(AECs)
  # 标准化数据
  AECs <- NormalizeData(AECs)
  
  # 识别高变基因
  AECs <- FindVariableFeatures(AECs)
  
  # 缩放数据
  AECs <- ScaleData(AECs)
  
  # PCA 降维
  AECs <- RunPCA(AECs)
  ElbowPlot(AECs)
  # 最近邻图和聚类
  AECs <- FindNeighbors(AECs, dims = 1:14)
  AECs <- FindClusters(AECs, resolution = 1)
  
  # UMAP 降维
  AECs <- RunUMAP(AECs, dims = 1:14)
  
  # 可视化新的聚类结果
  DimPlot(AECs, reduction = "umap", label = TRUE)
  
  AT_sub_markers <- FindAllMarkers(AECs, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  head(AT_sub_markers)
  head(AECs)
  
  # 生成颜色
  cell_type <- unique(AECs$seurat_clusters)
  colors <- colorRampPalette(brewer.pal(12, "Set3"))(length(cell_type))
  # 自定义配色
  DimPlot(AECs, reduction = "umap", group.by = "seurat_clusters", label = TRUE, label.size = 4) +
    labs(title = "UMAP OF AT2") +
    scale_color_manual(values = setNames(colors, cell_type))
  
  FeaturePlot(AECs, features = c("ESR1"))
  
  # 获取 Seurat clusters 的独特值
  condition <- unique(AECs$condition)
  colors <- colorRampPalette(brewer.pal(2, "Set3"))(length(condition))
  
  # 提取每个基因的单独 VlnPlot
  genes <- c("CEPT1")
  plots <- lapply(genes, function(gene) {
    VlnPlot(AECs, features = gene, group.by = "condition") +
      scale_fill_manual(values = colors) +
      labs(title = paste("Expression of", gene)) +
      theme(
        panel.background = element_blank(),
        panel.grid = element_blank(),
        title = element_text(size = 8, colour = "grey0"),
        axis.ticks = element_line(linewidth = 0.35, colour = "grey0"),
        text = element_text(size = 8, colour = "grey0"),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1, colour = "grey0"),
        axis.text.y = element_text(size = 8, colour = "grey0"),
        axis.title = element_blank(),
        legend.position = "none"
      )
  })
  
  # 使用 patchwork 拼接所有图形
  library(patchwork)
  final_plot <- wrap_plots(plots, ncol = 3)  # 设置列数
  print(final_plot)
  
  # # 差异基因分析
  at_markers <- FindMarkers(seurat, ident.1 = "AT2 cell", ident.2 = NULL, group.by = "cell_type")
  # 查看结果
  head(at_markers)
  table(AECs$condition)
  
  saveRDS(AECs, file = "AECs.rds")
  
  AECs <- readRDS("AECs.rds")
  if(T){
    # 提取AT2 注释
    # 初始化所有 cluster 的注释为原始编号
    cluster_anno <- as.character(Idents(AECs))
    
    # 手动定义部分 cluster 的细胞类型
    cluster_anno[Idents(seurat) %in% c(0,1,2,3,4,5,6,7,8,9,10,11,16,17)] <- "AT2"  # Cluster 0
    cluster_anno[Idents(seurat) %in% c(12,13,14,15)] <- "AT1"    # Cluster 1
    # cluster_annotations[Idents(seurat) == 2] <- "AT1 cell"    # Cluster 2
    # cluster_annotations[Idents(seurat) == 7] <- "AM"    # Cluster 7
    # 其他 cluster 暂不定义，保留原始编号
    
    # 分配细胞类型
    AECs$cell_type <- cluster_anno
    
    library(RColorBrewer)
    
    # 生成颜色
    cell_types <- unique(AECs$cell_type)
    colors <- c("skyblue","tomato")
    # 自定义配色
    DimPlot(AECs, reduction = "umap", group.by = "cell_type", label = TRUE, label.size = 4) +
      labs(title = "UMAP with Annotated Cell Types") +
      scale_color_manual(values = setNames(colors, cell_types))
    
    # 提取 AT2 cell
    seurat <- subset(AECs, subset = cell_type == "AT2")
    # 查看提取的结果
    print(seurat)
    # 标准化数据
    seurat <- NormalizeData(seurat)
    
    # 识别高变基因
    seurat <- FindVariableFeatures(seurat)
    
    # 缩放数据
    seurat <- ScaleData(seurat)
    
    # PCA 降维
    seurat <- RunPCA(seurat)
    ElbowPlot(seurat)
    # 最近邻图和聚类
    seurat <- FindNeighbors(seurat, dims = 1:15)
    seurat <- FindClusters(seurat, resolution = 0.5)
  }
  
  
  if(T){
    # 差异表达分析----
    DefaultAssay(AEC_obj) <- "RNA"
    
    # 添加分组信息----
    seurat$condition <- ifelse(seurat$phenotype == "ANA", "Control", "Asthma")
    seurat$condition <- factor(seurat$condition, levels = c("Control", "Asthma"))
    table(seurat$condition)
    # 设置分组
    Idents(AEC_obj) <- "condition"
    table(Idents(AEC_obj))
    # 进行差异表达分析
    AEC_obj_DEG <- FindMarkers(
      object = AEC_obj,
      # features = "CEPT1",
      ident.1 = "Asthma",        # 组1
      ident.2 = "Control",        # 组2
      min.pct = 0,             # 最小表达细胞比例
      test.use = "wilcox",        # 默认
      logfc.threshold = 0      # 最小 log2 fold-change 阈值
    )
    
    AEC_obj_DEG$gene <- rownames(AEC_obj_DEG) 
    write.csv(AEC_obj_DEG,file = "AEC_obj_DEG.csv", row.names = T)
    if(T){
      expression_data <- as.data.frame(FetchData(seurat_obj, 
                                                 vars = c("CEPT1","condition")))
      library(tidyr)
      long_predata <- expression_data %>% 
        pivot_longer(
        cols = !condition,      # 匹配以 "Gender_" 开头的列
        names_to = "name",             # 新的变量列名称
        values_to = "value"            # 新的值列名称
      )
      table(long_predata$name)
      # # 显示转换后的长数据格式的数据框
      print(long_predata)
      # 检查数据
      head(expression_data)
      
      for (name in unique(expression_data$name))
      {
        df2 <- long_predata[grep(paste0("^",name,"$"),long_predata$name, ignore.case = TRUE), ] 
        
        # 小提琴图
        long_predata <- long_predata %>% mutate(condition = factor(condition, levels = c("Control", "Asthma")))
        long_predata <- arrange(long_predata, condition)
        
        # 移除表达为0的样本
        # data <- data[data$HTRA2 > 0, ]
        # 定义一个函数来去除每组的异常值
        # if(T){
        #   # 定义去除异常值的函数，允许设置自定义的 IQR 倍数
        #   remove_outliers <- function(x, multiplier = 1) {
        #     Q1 <- quantile(x, 0.25)
        #     Q3 <- quantile(x, 0.75)
        #     IQR <- Q3 - Q1
        #     # 使用自定义的 IQR 倍数返回是否为非异常值的逻辑向量
        #     x >= (Q1 - multiplier * IQR) & x <= (Q3 + multiplier * IQR)
        #   }
        #   
        #   # 使用 dplyr 包按组过滤异常值，并在结果中标记是否为异常值
        #   library(dplyr)
        #   data_clean <- data %>%
        #     group_by(group) %>%
        #     mutate(is_outlier = !remove_outliers(value, multiplier = 1)) %>% # 标记异常值
        #     filter(!is_outlier) %>% # 仅保留非异常值
        #     ungroup()
        #   
        #   data <- data_clean
        #  }
        # 绘制箱线图----
        table(long_predata$condition)
        p <- ggplot(long_predata, aes(condition, value, fill = condition)) +
          geom_violin(trim = FALSE, alpha = 0.3, scale = "width", linewidth = 0.35) + # 增加小提琴图，透明度调低
          geom_boxplot(aes(fill = condition),
                       alpha = 0.5,
                       linewidth = 0.5,
                       notch = FALSE,
                       width = 0.5,
                       outlier.size = 0.05,
                       outlier.color = "black"
          ) +
          geom_jitter(mapping = aes(colour = condition),
                      shape = 16,
                      position = position_jitter(0.2),
                      size = 0.3,
                      alpha = 0.8
          ) +
          geom_signif(
            comparisons = list(c("Control", "Asthma")), 
            map_signif_level = T, 
            test = "t.test", # 或 "wilcox.test" 视具体情况
            step_increase = 1,
            tip_length = 0,
            textsize = 3.5,
            size = 0.35
          ) +
          stat_summary(fun = mean, geom = "point", color = "blue", size = 2) +
          stat_summary(fun = median, geom = "point", color = "red", size = 2) +
          labs(
            title = name,
            x = "Condition",
            y = bquote(italic(.(name)) ~ " Relative Expression")
          ) +
          theme(
            # panel.grid.major = element_line(color = "grey90", size = 0.35),
            # panel.border = element_rect(color = "black", fill = NA, size = 0.35),
            panel.background = element_blank(),
            panel.grid = element_blank(),
            axis.line = element_line(linewidth = 0.35),
            axis.ticks = element_line(linewidth = 0.35, colour = "grey0"),
            text = element_text(size = 8, colour = "grey0"),
            axis.text.x = element_text(size = 8,  colour = "grey0"),
            axis.text.y = element_text(size = 8, colour = "grey0"),
            axis.title.y = element_text(size = 8, colour = "grey0"),
            axis.title.x = element_blank(),
            plot.title = element_text(size = 8, face = "bold"),
            legend.position = "none"
          ) +
          scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
          # scale_y_log10(expand = expansion(mult = c(0, 0.1))) +
          scale_fill_manual(values = c("Control" = "#87cEEB", "Asthma" = "#FF5722")) +
          scale_colour_manual(values = c("Control" = "#87cEEB", "Asthma" = "#FF5722"))
       
        # # 创建 label_data，按 cell_type 和 condition 分组计算表达比例
        # label_data <- long_predata %>%
        #   group_by(condition) %>%
        #   summarise(rate = mean(name > 0) * 100, .groups = "drop")
        # 
        # # 添加文字标签到图下方（y = -0.5）
        # p <- p + 
        #   geom_text(
        #     data = label_data,
        #     aes(x = condition, y = -0.5, label = paste0(round(rate, 1), "%"), group = condition, color = condition),
        #     position = position_dodge(width = 0.8),
        #     size = 3,
        #     show.legend = FALSE
        #   )
        # 显示图
        print(p)
        ggsave(paste0("HTRA2", ".tiff"), plot = p, height = 1.8, width = 2,dpi = 300)
        ggsave(paste0("HTRA2", ".pdf"), plot = p, height = 1.8, width = 2,dpi = 300)
        # 箱线图，一起绘制----
        p <- ggplot(long_predata, aes(name, value, fill = condition)) + 
          geom_violin(trim = FALSE, alpha = 0.25, scale = "width", linewidth = 0.35, 
                      position = position_dodge(width = 0.88)) +  # 确保不同组分开
          geom_boxplot(
            alpha = 0.3,  
            linewidth = 0.45,
            notch = FALSE,
            width = 0.6,
            outlier.shape = NA,
            position = position_dodge(width = 0.88)  # 使箱线图对齐
          ) +
          geom_jitter(
            aes(colour = condition),  # 统一颜色映射，避免生成额外 legend
            shape = 16,
            position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.88),  # 避免数据混在一起
            size = 0.5,
            alpha = 0.6
          ) +
          theme_bw() + 
          labs(x = NA, y = "Log-normalized Expression") +
          theme(
            legend.position = "right",  
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 0, hjust = 0.5, size = 8, fAEC = "italic", colour = 'black'),
            axis.text.y = element_text(size = 8, colour = 'black')
          ) +
          scale_fill_manual(
            values = c("Control" = "royalblue", "Asthma" = "tomato"),
            breaks = c("Control", "Asthma"),
            labels = c("Control" = "Control", "Asthma" = "Asthma")
          ) +
          scale_color_manual(  # 统一 jitter 颜色
            values = c("Control" = "royalblue", "Asthma" = "tomato")
            ) +
          stat_compare_means(
            aes(group = condition),
            label = "p.format",
            size = 3,
            method = "wilcox.test" #t.test
          ) 
        p  
        # 保存图片,多个柱子，宽度加0.5
        ggsave(paste0("核心基因-单细胞", ".tiff"), plot = p, height = 3, width = 10,dpi = 300)
        ggsave(paste0("核心基因-单细胞", ".pdf"), plot = p, height = 3, width = 10,dpi = 300)
        
        result <- t.test(value~condition, data, var.equal = T)#方差齐
        print(result)
        
      }
    }
    
    # 单基因相关性----
    # 提取 CDKN1A 和 YAP1 的表达量
    {
      seurat$CDKN1A_expr <- FetchData(seurat, vars = "CDKN1A")
      seurat$YAP1_expr <- FetchData(seurat, vars = "YAP1")
      
      # 绘制散点图
      library(ggplot2)
      # 过滤低表达的细胞
      filtered_cells <- subset(seurat, subset = CDKN1A > 0 | YAP1 > 0)
      p6 <- ggplot(filtered_cells@meta.data, aes(x = CDKN1A_expr, y = YAP1_expr, color = condition)) +
        geom_point(aes(color = condition), alpha = 0.6, size = 2) +  
        scale_color_manual(values = c("Healthy" = "skyblue", "Disease" = "tomato")) + 
        labs(
          title = "Correlation between YAP1 and ABLIM1",
          x = "CDKN1A Expression",
          y = "YAP1 Expression",
          color = "Condition"
        ) +
        stat_smooth(color = "grey30", formula = y ~ x, fill = "lightgrey", method = "lm") +
        stat_cor(
          # aes(color = condition),
          method = "pearson",  # 或替换为 "spearman"
          digits = 2  # 保留两位小数
        ) +
        theme(
          # panel.grid.major = element_line(color = "grey90", size = 0.35),
          # panel.border = element_rect(color = "black", fill = NA, size = 0.35),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(linewidth = 0.35),
          axis.ticks = element_line(linewidth = 0.35, colour = "grey0"),
          text = element_text(size = 12, colour = "grey0"),
          axis.text.x = element_text(size = 8,  colour = "grey0"),
          axis.text.y = element_text(size = 8, colour = "grey0"),
          axis.title.y = element_text(size = 8, colour = "grey0"),
          axis.title.x = element_text(size = 8, colour = "grey0"),
          plot.title = element_text(size = 8, fAEC = "bold"),
          legend.position = "right"
        )
      p6
      
      ggsave(paste0("YAP1 and ABLIM1", ".tiff"), plot = p6, height = 3, width = 4,dpi = 300)
      ggsave(paste0("YAP1 and ABLIM1", ".pdf"), plot = p6, height = 3, width = 4)
    }
    
    
    library(dplyr)
    # 添加表达分类标签
    filtered_cells@meta.data <- filtered_cells@meta.data %>%
      mutate(
        Expression_Category = case_when(
          ITGAX_expr > 0 & PIEZO1_expr > 0 ~ "Both Expressed",
          ITGAX_expr > 0 & PIEZO1_expr == 0 ~ "CDKN1A Only",
          ITGAX_expr == 0 & PIEZO1_expr > 0 ~ "YAP1 Only",
          TRUE ~ "None"
        )
      )
    
    # 查看分类的细胞数量
    table(filtered_cells@meta.data$Expression_Category)
    
    # 绘制堆叠条形图
    ggplot(filtered_cells@meta.data, aes(x = Expression_Category, fill = condition)) +
      geom_bar(position = "fill") +
      labs(
        title = "Expression Patterns of YAP1 and CDKN1A",
        x = "Expression Category",
        y = "Proportion",
        fill = "Condition"
      ) +
      theme_minimal()
    
    
    ## 多基因相关性----
    if(T){
      
      # 查看当前 Seurat 对象的默认 Assay
      DefaultAssay(seurat) <- "RNA"
      
      
      # 查看所有 Assay
      Assays(seurat)
      
      # 查看 RNA Assay 的所有层
      Layers(seurat[["RNA"]])
      
      target_genes <- c("CDKN1A", "CDKN2A", "YAP1", "ABLIM1", "TEAD1", "TEAD2", "TEAD3", "TEAD4")
      missing_genes <- setdiff(target_genes, rownames(seurat))
      if (length(missing_genes) > 0) {
        stop("以下基因不存在于数据中：", paste(missing_genes, collapse = ", "))
      } else {
        cat("所有基因都存在于数据中。\n")
      }
      
      # 获取 RNA 层的表达矩阵
      expr_matrix <- GetAssayData(seurat, assay = "RNA", slot = "data")
      
      # 提取目标基因
      filtered_genes <- expr_matrix[target_genes, ]
      
      # 查看过滤后的基因表达矩阵
      filtered_genes <- filtered_genes[rowSums(filtered_genes) > 0, ]
      
      # 创建包含目标基因的 Seurat 对象
      filtered_cells <- CreateSeuratObject(filtered_genes)
      
      # 计算相关性矩阵
      cor_matrix <- cor(as.matrix(filtered_genes), method = "pearson")
      
      # 绘制热图
      library(pheatmap)
      pheatmap(
        cor_matrix,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        display_numbers = TRUE,
        color = colorRampPalette(c("blue", "white", "red"))(100)
      )
      # 确保基因存在
      target_genes <- c("YAP1", "ABLIM1", "TEAD1", "TEAD2", "TEAD3", "TEAD4", "CDKN1A", "CDKN2A")
      # 提取 RNA 数据的归一化表达矩阵
      expr_matrix <- GetAssayData(seurat, assay = "RNA", slot = "data")
      
      # 过滤目标基因
      filtered_genes <- expr_matrix[target_genes, ]
      
      # 使用基因表达的逻辑条件筛选细胞
      filtered_cells <- subset(seurat, subset = YAP1 > 0 | ABLIM1 > 0 | TEAD1 > 0 | TEAD2 > 0 | TEAD3 > 0 | TEAD4 > 0 | CDKN1A > 0 | CDKN2A > 0)
      # 提取目标基因的表达矩阵
      expr_matrix_filtered <- GetAssayData(filtered_cells, assay = "RNA", slot = "data")[target_genes, ]
      # # 转换稀疏矩阵为普通矩阵
      # expr_matrix_filtered_dense <- as.matrix(expr_matrix_filtered)
      # 定义两组基因
      group1 <- c("YAP1", "ABLIM1", "TEAD1", "TEAD2", "TEAD3", "TEAD4")
      group2 <- c("CDKN1A", "CDKN2A")
      
      # 确保基因名存在
      group1 <- group1[group1 %in% rownames(expr_matrix_filtered)]
      group2 <- group2[group2 %in% rownames(expr_matrix_filtered)]
      
      # 提取两组基因的表达矩阵
      expr_group1 <- as.matrix(expr_matrix_filtered[group1, ])  # 第一组基因
      expr_group2 <- as.matrix(expr_matrix_filtered[group2, ])  # 第二组基因
      
      # 计算相关性：行是基因，列是细胞，因此不需要转置
      cor_matrix <- cor(t(expr_group1), t(expr_group2), method = "pearson")
      # 转置矩阵以便热图展示时第一组基因为行，第二组基因为列
      cor_matrix <- t(cor_matrix)
      # 查看相关性矩阵
      print(cor_matrix)
      # # 定义显著性阈值
      # threshold <- 1
      # 
      # # 创建显著性矩阵：相关性大于阈值标记为 TRUE
      # significance_matrix <- abs(cor_matrix) > threshold
      # 
      # # 替换显著性矩阵的值为星号或空字符串
      # significance_labels <- ifelse(significance_matrix, "*", "")
      # 
      # 定义颜色渐变
      # color_palette <- colorRampPalette(c("lightblue", "white", "tomato"))(100)
      
      # 绘制热图
      pheatmap(cor_matrix, 
               main = "Correlation: Group1 vs Group2", 
               cellwidth = 30, cellheight = 50, fontsize = 10,  # 控制格子大小和字体
               # display_numbers = significance_labels,  # 显示显著性标记
               number_color = "black",  # 星号颜色
               cluster_rows = FALSE, cluster_cols = FALSE)  # 设置颜色渐变
      
      
    }
    # 查看结果
    head(at2_diff_gene)
    at2_diff_gene$gene <- row.names(at2_diff_gene)
    # 热图，差异可视化
    top_genes <- rownames(head(at2_diff_gene, 100))  # 取前20个显著基因
    colors <- c("Healthy" = "skyblue", "Disease" = "tomato")
    # 确保 "condition" 是一个因子，并设置因子水平
    # 检查 scale.data 中是否包含所有基因
    scale_genes <- rownames(GetAssayData(seurat, layer = "scale.data"))
    missing_genes <- setdiff(top_genes, scale_genes)
    
    # 如果有缺失基因，对它们进行标准化
    if (length(missing_genes) > 0) {
      seurat <- ScaleData(seurat, features = union(scale_genes, top_genes))
    }
    
    # 确保最终的基因列表存在于 scale.data 中
    top_genes <- intersect(top_genes, rownames(GetAssayData(seurat, layer = "scale.data")))
    
    p5 <- DoHeatmap(seurat, 
                    features = top_genes, 
                    group.by = "condition",
                    slot = "scale.data",
                    label = F,
                    group.colors = colors,
                    # disp.min = -2.5,       # 热图最小值
                    # disp.max = 2.5,         # 热图最大值
                    size = 6,       # 调整基因名称字体大小
                    angle = 0       ) + # 列标签旋转角度
      labs(title = "Top Differentially Expressed Genes")
    p5
    ggsave(paste0("AT2_hm", ".tiff"), plot = p5, height = 10, width = 8,dpi = 300)
    ggsave(paste0("AT2_hm", ".pdf"), plot = p5, height = 10, width = 8)
    # 提取 scale.data 中的基因名称
    scale_genes <- rownames(GetAssayData(seurat, layer = "scale.data"))
    
    # 检查缺失的基因
    missing_genes <- setdiff(top_genes, scale_genes)
    missing_genes
    
    
    # 差异表达分析
    markers_select <- FindMarkers(seurat, ident.1 = "Healthy", ident.2 = "Disease", features = "ABLIM1")
    print(markers)
    # 调整因子顺序，健康在左，疾病在右
    seurat$condition <- factor(seurat$condition, levels = c("Healthy", "Disease"))
    
    VlnPlot(seurat, features = c("ABLIM1","YAP1","TEAD1","TEAD2"), group.by = "condition") #初步查看
    
    data <- as.data.frame(FetchData(seurat, vars = c("ABLIM1","YAP1","TEAD1","TEAD2","TEAD3","TEAD4", "condition")))
    data <- as.data.frame(FetchData(seurat, vars = c("CDKN1A","CDKN2A","TP53","SFTPA1", "condition")))
    
    # 检查数据
    head(data)
    
    library(tidyr)
    long_predata <- pivot_longer(
      data = data,
      cols = !c(condition),      # 匹配以 "Gender_" 开头的列
      names_to = "name",             # 新的变量列名称
      values_to = "value"            # 新的值列名称
    )
    table(long_predata$name)
    # # 显示转换后的长数据格式的数据框
    print(long_predata)
    
    for (gene in unique(long_predata$name)) {
      data <- long_predata[grep(paste0("^",gene,"$"),long_predata$name, ignore.case = TRUE), ] 
      data <- data %>% mutate(group = factor(condition, levels = c("Healthy", "Disease")))
      data <- arrange(data, group)
      # 移除表达为0的样本
      # data <- data[data$value > 0, ]
      # 绘制箱线图
      p <- ggplot(data, aes(x = group, y = value, fill = group)) +
        geom_violin(trim = FALSE, alpha = 0.3, scale = "width", linewidth = 0.35) + # 增加小提琴图，透明度调低
        geom_boxplot(aes(fill = group),
                     alpha = 0.3,
                     linewidth = 0.35,
                     notch = FALSE,
                     width = 0.5,
                     outlier.size = 0.05,
                     outlier.color = "black"
        ) +
        geom_jitter(mapping = aes(colour = group),
                    shape = 16,
                    position = position_jitter(0.2),
                    size = 0.8,
                    alpha = 0.6
        ) +
        geom_signif(
          comparisons = list(c("Disease", "Healthy")), 
          map_signif_level = T, 
          test = "t.test", # 或 "wilcox.test" 视具体情况
          step_increase = 1,
          tip_length = 0,
          textsize = 3.5,
          size = 0.35
        ) +
        stat_summary(fun = mean, geom = "point", color = "blue", size = 2) +
        stat_summary(fun = median, geom = "point", color = "red", size = 2) +
        labs(
          title = gene, 
          x = "Group", 
          y = bquote(italic(.(gene)) ~ "Relative Expression")
        ) +
        theme_minimal() +
        theme(
          # panel.grid.major = element_line(color = "grey90", size = 0.35),
          # panel.border = element_rect(color = "black", fill = NA, size = 0.35),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(linewidth = 0.35),
          axis.ticks = element_line(linewidth = 0.35, colour = "grey0"),
          text = element_text(size = 8, colour = "grey0"),
          axis.text.x = element_text(size = 8,  colour = "grey0"),
          axis.text.y = element_text(size = 8, colour = "grey0"),
          axis.title.y = element_text(size = 8, colour = "grey0"),
          axis.title.x = element_blank(),
          plot.title = element_blank(),
          legend.position = "none"
        ) +
        scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
        # scale_y_log10(expand = expansion(mult = c(0, 0.1))) +
        scale_fill_manual(values = c("Healthy" = "#87cEEB", "Disease" = "#FF5722")) +
        scale_colour_manual(values = c("Healthy" = "#87cEEB", "Disease" = "#FF5722"))
      # 显示图
      print(p)
      # 保存图片,多个柱子，宽度加0.5
      ggsave(paste0(gene, ".tiff"), plot = p, height = 3, width = 3,dpi = 300)
      ggsave(paste0(gene, ".pdf"), plot = p, height = 3, width = 3,dpi = 300)
    }
  }
  
  saveRDS(seurat, file = "seurat.rds")
  save(seurat, at2_markers,file = "AT2_select.Rdata")
  table(seurat$cell_type)
  table(seurat$sample, seurat$cell_type)
  
  
  
  result <- t.test(value~group, data, var.equal = T)#方差齐
  print(result)
}

# 富集分析----

load("AT2_select.Rdata")

# 筛选差异基因
at2_diff_gene$gene <- rownames(at2_diff_gene)
sig_genes <- subset(at2_diff_gene, p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)
# 提取显著基因
gene_list <- sig_genes$gene

library(clusterProfiler)
library(org.Hs.eg.db)

# 差异基因转换为 ENTREZ ID
gene_entrez <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# EnrichGO富集----
##CC表示细胞组分，MF表示分子功能，BP表示生物学过程，ALL表示同时富集三种过程，选自己需要的,我一般是做BP,MF,CC这3组再合并成一个数据框，方便后续摘取部分通路绘图。
enrichGO_ALL <- enrichGO(gene = gene_entrez$ENTREZID,#我们上面定义了
                         OrgDb=org.Hs.eg.db,    #org.Hs.eg.db, org.Mm.eg.db
                         keyType = "ENTREZID",
                         ont = "ALL",#富集的GO类型
                         pAdjustMethod = "BH",#这个不用管，一般都用的BH
                         minGSSize = 1,
                         pvalueCutoff = 0.05,#P值可以取0.05
                         qvalueCutoff = 0.3,
                         readable = TRUE #TURE表示以基因名显示
) 

# dim(enrichGO_ALL[enrichGO_ALL$ONTOLOGY=='BP',])
# dim(enrichGO_ALL[enrichGO_ALL$ONTOLOGY=='CC',])
# dim(enrichGO_ALL[enrichGO_ALL$ONTOLOGY=='MF',])
dim(enrichGO_ALL)

# 去除GO富集term中高度相似的冗余
go.filter <- clusterProfiler::simplify(enrichGO_ALL,cutoff = 0.5,
                                       by = "p.adjust",
                                       select_fun = min)
dim(go.filter)
go.filter_res <- go.filter@result
go.filter_res <- go.filter_res[order(go.filter_res$RichFactor, decreasing = TRUE),]
write.csv(go.filter_res,file = "EnrichGO.filter.csv",row.names = T)

barplot(go.filter,drop = TRUE,
        showCategory = 6,split = "ONTOLOGY") +
  fAECt_grid(ONTOLOGY~., scale ='free')
# go.filter_res <- mutate(go.filter_res, RichFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))#添加RichFactor列

#分别提取BP,CC,MF数据
result_BP <- go.filter_res %>% filter(ONTOLOGY=='BP')
result_CC <- go.filter_res %>% filter(ONTOLOGY=='CC')
result_MF <- go.filter_res %>% filter(ONTOLOGY=='MF')
#取前10行
BP<-result_BP[1:10,]
CC<-result_CC[1:10,]
MF<-result_MF[1:10,]
all <- rbind(BP,CC,MF)


##开始绘制GO柱状图
# 拆分分数形式的字符串，计算分子除以分母，sapply函数的特点是一次性对一堆数据执行某个函数2，并且返回值为向量或矩阵
all$GeneRatio <- sapply(strsplit(all$GeneRatio, "/"), 
                        function(x) as.numeric(x[1]) / as.numeric(x[2]))
# 横着的柱状图
if(T){
  # color <- c("#66C3A5", "#8DA1CB", "#FD8D62") # 定义ONTOLOGY对应的颜色
  # ontology_colors <- c("BP" = "#66C3A5", "CC" = "#8DA1CB", "MF" = "#FD8D62")
  pic <- ggplot(all, aes(x = RichFactor, y = reorder(Description, RichFactor))) + #横纵轴取值
    geom_point(aes(size=Count,fill=-log(p.adjust)),
               color = "black",
               shape=21) + #柱状图的宽度，可以自己设置
    
    scale_fill_gradient(low = "#4169E1", high = "tomato") +
    scale_color_manual(values = color) + ###颜色
    fAECt_grid(ONTOLOGY~.,
               scale='free_y',spAEC = 'free_y'
    )+
    xlab("RichFactor") +
    ylab("GO term") +
    labs(title = "The Most Enriched GO Terms",)+
    theme_bw()+
    theme(panel.grid =element_blank(),
          panel.border = element_rect(linewidth = 0.35),
          panel.grid.major = element_line(color = "gray", linewidth = 0.35, linetype = "dashed"),
          plot.title = element_text(size = 8,hjust=0.5),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 8),
          text = element_text(family = "sans", size = 12),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_text(size = 8),
          axis.text = element_text(size = 8,colour = "grey0"),
          strip.text = element_text(size = 8, colour = "grey0"),
          strip.background = element_rect(fill = "skyblue",linewidth = 0.35)
    )
  
  
  pic
  ggsave("GO.tiff", plot = pic, height = 10, width = 8,dpi = 300)
  ggsave("GO.pdf", plot = pic, height = 10, width = 8,dpi = 300)
}

# EnrichKEGG----
if(T) {
  rich_kk <- enrichKEGG(gene = gene_entrez$ENTREZID,
                        organism = "hsa",
                        keyType = "kegg", #kEGG数据库
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05)
  rich_kk_res <- setReadable(rich_kk, 'org.Hs.eg.db', 'ENTREZID')
  rich_kk_res <- rich_kk_res@result#自己记得保存结果哈！
  write.csv(rich_kk_res, file = "enrich_kk_results.csv", row.names = T)
  # # 自定义排序
  # rownames(rich_kk_res) <- 1:nrow(rich_kk_res)#行名设置为了从 1 到行数的整数
  # rich_kk_res$description <- factor(rev(as.integer(rownames(rich_kk_res))),
  #                                   labels = rev(sub("(.+)-.*", "\\1", 
  #                                                    rich_kk_res$Description)))
  # 按GeneRatio排序
  rich_kk_res <- mutate(rich_kk_res, geneRatio = parse_ratio(GeneRatio)) %>%
    arrange(desc(geneRatio))
  # 按p.adjust升序
  rich_kk_res <- rich_kk_res[order(rich_kk_res$p.adjust), 
                             ,drop = FALSE]
  
  rich_kk_res <- rich_kk_res[order(rich_kk_res$RichFactor, decreasing = TRUE), 
                             ,drop = FALSE]
  dotplot(rich_kk)
  # KEGG, Dotplot 自定义气泡图
  library(tidyverse)
  # 删除不要的命名部分
  # rich_kk_res$Description <- sub(" - Mus musculus \\(house mouse\\)$", "", rich_kk_res$Description)
  # rich_kk_res$GeneRatio <- as.numeric(str_split(rich_kk_res$GeneRatio, pattern = '/', simplify = T)[,1])
  # rich_kk_res <- mutate(rich_kk_res, RichFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))#添加RichFactor列
  
  if(T){
    pic <- ggplot(rich_kk_res[1:20,], aes(x = RichFactor,y = fct_reorder(Description, RichFactor))) + #横纵轴取值
      geom_point(aes(size=Count,fill=p.adjust),
                 color = "black",
                 shape=21) + #柱状图的宽度，可以自己设置
      # scale_fill_manual(values = color) + ###颜色
      scale_fill_gradient(low = "#FF6347", high = "#4169E1") +
      # coord_flip() + ##这一步是让柱状图横过来，不加的话柱状图是竖着的
      ylab("KEGG term") +
      xlab("RichFactor") +
      labs(title = "The Most Enriched KEGG Terms")+
      # fAECt_grid(rows = vars(subcategory),scales = 'free_y',spAEC = 'free_y') +
      theme_bw()+
      theme(panel.grid =element_blank(),
            panel.border = element_rect(linewidth = 0.35),
            panel.grid.major = element_line(color = "gray", linewidth = 0.35, linetype = "dashed"),
            plot.title = element_text(size = 8,hjust=0.5),
            legend.title = element_text(size = 8),
            legend.text = element_text(size = 8),
            legend.background = element_rect(fill = 'transparent'),
            text = element_text(family = "sans", size = 12),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_text(size = 8),
            axis.text = element_text(size = 8,colour = "grey0")
      )#angle是坐标轴字体倾斜的角度，可以自己设置
    pic
    ggsave("KEGG_category.tiff", plot = pic, height = 10, width = 8,dpi = 300)
    ggsave("KEGG_category.pdf", plot = pic, height = 10, width = 8,dpi = 300)
  }
  dev.off()
  
  # 分亚组----
  pic <- ggplot(rich_kk_res[1:20,], aes(x = RichFactor, y = fct_reorder(Description, RichFactor),fill = category)) +
    geom_bar(stat = "identity", width = 0.8) + #柱状图的宽度，可以自己设置
    scale_fill_manual(values = c("#FF6347","#B3DE69","#E6AB02","#80B1D3","#BEBADA")) +
    
    # scale_fill_manual(values = c("#FF6347","#B3DE69","#BEBADA")) +
    # fAECt_grid(rows = vars(category),scales = 'free_y',spAEC = 'free_y') +
    # ylab("KEGG term") +
    xlab("RichFactor") +
    labs(title = "The Most Enriched KEGG Terms",
         # x = bquote(~-Log[10]~italic("P-value")),
         # x = "GeneRatio",
         y = "KEGG term",
         fill  = "Category") + #离散型用col
    geom_text(aes(label = GeneRatio),size = 3) +
    # geom_hline(yintercept = -log10(as.numeric(0.05)),color = 'grey',size = 0.35,lty = 'dashed') + #添加虚线
    # geom_vline(xintercept = -log10(as.numeric(0.3)),color = 'grey',size = 0.35,lty = 'dashed') +
    # scale_x_continuous(expand = c(0,0),breaks = seq(0,0.2,5), limits = c(0,0.2,5)) + 
    # scale_y_discrete(labels = function(dat) str_wrap(rich_kk_res[1:30,],width = 25)) +
    theme_bw() + theme(panel.grid =element_blank(),
                       panel.border = element_rect(linewidth = 0.35),
                       # panel.grid.major = element_line(color = "gray", linewidth = 0.35, linetype = "dashed"),
                       plot.title = element_text(size = 8,hjust=0.5),
                       legend.title = element_text(size = 8),
                       legend.text = element_text(size = 8),
                       legend.background = element_rect(fill = 'transparent'),
                       text = element_text(family = "sans", size = 12),
                       axis.line = element_blank(),
                       axis.ticks = element_line(linewidth = 0.35,colour = "grey0"),
                       axis.title = element_text(size = 8),
                       axis.text = element_text(size = 8,colour = "grey0")
    ) 
  pic
  ggsave("KEGG_category_1.tiff", plot = pic, height = 10, width = 8,dpi = 300)
  ggsave("KEGG_category_1.pdf", plot = pic, height = 10, width = 8,dpi = 300)
  # 映射通路----
  ###设置工作目录到想要的地方
  # dir.create("C:/wj/RNASeq/pathwayview_out")
  setwd("C:/wj/RNASeq/pathwayview_out")
  
  select_pathway <- c("mmu04020","mmu04060","mmu04151","mmu04062","mmu04514","mmu04630","mmu04657","mmu04613","mmu05418") #选择所需通路的ID号
  pathview(gene.data     = geneList,
           pathway.id    = select_pathway,
           out.suffix = "map",
           species       = 'mmu' ,      # 人类hsa 小鼠mmu 
           kegg.native   = T, #TRUE输出完整pathway的png文件，F输出基因列表的pdf文件 
           new.signature = F, #pdf是否显示pathway标注
           limit         = list(gene=2.5, cpd=1)
  )
  
}

if(T){#GSEA_GO####
  gseGO <- gseGO(geneList = geneList,
                 OrgDb = org.Mm.eg.db, 
                 ont = "all", #可选条目BP/CC/MF
                 pAdjustMethod = "BH", #p值的校正方式
                 pvalueCutoff = 0.05) # pvalue的阀值) #是否将entrez id转换为symbol
  tmp <- setReadable(gseGO, 'org.Mm.eg.db', 'ENTREZID')
  gseGO_res <- tmp@result
  gseGO_res <- gseGO_res %>% arrange(desc(enrichmentScore))
  write.csv(gseGO_res,file = "gseGO.csv",row.names = T)
  
  # 简要查看----
  dotplot(go_res)
  dotplot(
    gseGO,
    color = "p.adjust",
    showCategory=10,
    split=".sign") + fAECt_grid(.~.sign)
  
  # GSEA_GO, Dotplot 自定义气泡图----
  library(tidyverse)
  Dotplot <- function(gseGO_res, title="") {
    dotplot <- ggplot(cbind(gseGO_res, Order = nrow(gseGO_res):1)) +
      geom_point(mapping = aes(x = -log10(p.adjust), y = Order, 
                               size = setSize, fill = enrichmentScore),
                 shape = 21) + 
      scale_fill_gradientn(colours = c("blue", "red")) + #自定义配色
      scale_y_continuous(position = "left", 
                         breaks = 1:nrow(gseGO_res), 
                         labels = Hmisc::capitalize(rev(gseGO_res$Description))) +
      # scale_x_continuous(breaks = c(3, 4,5,6),
      #                   #breaks = seq(0, xmax+5, 5),
      #                   limits = c(3,6),
      #                   expand = expansion(mult = c(.05, .05))) + #两边留空
      labs(x = "-Log10(p.adjust)", y = NULL) +
      guides(size = guide_legend(title = "setSize"),
             fill = guide_colorbar(title = "enrichmentScore")) +
      ggtitle(title) +
      theme_bw() +
      theme(panel.grid =element_blank(),
            panel.grid.major = element_line(color = "gray", linetype = "dashed"),
            panel.border = element_rect(color = "black", linewidth = 1),
            axis.text = element_text(size = rel(1)),
            title = element_text(size = rel(1))) #去除网格线
    # dotplot %>% ggplotGrob()#  %>% cowplot::plot_grid()
    # # fname=paste0(resultdir, '/', filemark,'.pdf')
    # # ggsave(fname, width = 8, height = 10)
    return(dotplot)
  }
  range(round(-log10(gseGO_res$p.adjust)))
  p_df <- gseGO_res[1:20,]
  Dotplot(p_df, title="GSEA-GO Enrich")
  ggsave(file = "Dotplot_GSEA-GO_top20_enrichmentscore.pdf", width = 8, height = 7)
  #gh <- ggplotGrob(dotplot)
  #gd <- ggplotGrob(dotplotk)
  #cowplot::plot_grid(gh, gd, rel_widths = c(1.2, 1))
  # 有向无环图
  plotGOgraph(go_res)#只能单一条目BP/CC/MF
  # 基因网络图
  cnetplot(gseGO_res,foldchange = geneList,
           circular = TRUE,
           color.params = list(edge = TRUE,
                               category = "#E5C494",
                               gene = "#B3B3B3"),
           showCategory = 10,
           cex.params = list(gene_label = 0.8,
                             category_label = 0.8,
                             category_node = 0.8,
                             gene_node = 0.8)
  )
  ggsave("go_res_circle.pdf",width=13,height=9)
  dev.off()
  ?cnetplot
}

# 转换基因ID
Genes <- bitr(at2_diff_gene[,6],
              fromType = "SYMBOL",#输入数据的类型
              toType = c("ENTREZID"),#要转换的数据类型
              OrgDb = org.Hs.eg.db)#物种
# 合并####
at2_diff_gene <-  merge(at2_diff_gene,Genes,
                        by.x="gene",
                        by.y="SYMBOL")
table(colnames(at2_diff_gene))
# 创建GSEA分析基因列表
GSEA_list <- at2_diff_gene %>% dplyr::select(ENTREZID,avg_log2FC)
# 降序
GSEA_lmn <- GSEA_list %>% dplyr::arrange(desc(avg_log2FC))#方法一
# GSEA_list = GSEA_list[order(GSEA_list$log2FoldChange,decreasing = TRUE),]#方法二
# GSEA_list
GSEA_lmn
geneList = GSEA_lmn[,2]
names(geneList) = as.character(GSEA_lmn[,1])
geneList
# GSEA_Kegg----
if(T)
{
  gseKK <- gseKEGG(geneList = geneList,
                   organism = "hsa",
                   pAdjustMethod = "BH",
                   minGSSize= 2,
                   maxGSSize= 500,
                   pvalueCutoff = 0.9,
                   verbose= FALSE)
  gseKK@readable
  tmp <- setReadable(gseKK,'org.Hs.eg.db', 'ENTREZID')
  gseKK_res <- tmp@result
  gseKK_res <- gseKK_res %>% arrange(desc(enrichmentScore))
  write.csv(gseKK_res, file = "gesKK.csv", row.names = T)
  
  tmp@result$Description <- sub(" - Mus musculus \\(house mouse\\)$", "", tmp@result$Description)
  
  # 排序
  gsea_res_order <- tmp[order(tmp$enrichmentScore,decreasing = TRUE),]
  # 绘图代码,GSEA 图
  library(enrichplot)
  if(T){
    geneSetIDs <- c("hsa04115")
    p <- gseaplot2(tmp, geneSetIDs,
                   title = "",
                   base_size = 10,
                   pvalue_table = T,
                   rel_heights = c(1.5,0.5,1),
                   color = colorspAEC::rainbow_hcl(3),
                   ES_geom = "line")
    
    p[[1]] <- p[[1]] +  
      # scale_color_viridis_d(labels = tmp$Description)+ #颜色参数
      geom_hline(yintercept = 0,
                 color = "grey75", 
                 linewidth = 0.35, 
                 linetype = 2)
    theme(legend.direction = "horizontal",
          text = element_text(size = 12))
    # p[[2]] <- p[[2]] +
    # scale_color_viridis_d(labels = tmp$Description)#颜色参数
    p[[3]] <- p[[3]] +
      geom_hline(yintercept = 0, color = "steelblue", 
                 linewidth = 0.35,
                 linetype = 2)
    p
    ggsave("GSEA_KEGG_p53.pdf",width=8,height=6)
    ggsave("GSEA_KEGG_p53.tiff",width=8,height=6)
  }
  
  # 通路活性分析----
  # 提取通路基因集
  pathway1_genes <- unlist(strsplit(gseKK_res[rownames(gseKK_res) == "hsa04218", "core_enrichment"], "/")) #p53
  pathway2_genes <- unlist(strsplit(gseKK_res[rownames(gseKK_res) == "hsa04390", "core_enrichment"], "/")) #Hippo signaling pathway
  
  
  library(GSVA)
  
  # 定义基因集列表
  pathway_list <- list(
    Pathway1 = pathway1_genes,
    Pathway2 = pathway2_genes
  )
  
  # 提取表达矩阵
  expr_matrix <- GetAssayData(seurat, assay = "RNA", slot = "data")
  
  # 筛选基因集
  pathway_list <- lapply(pathway_list, function(genes) {
    genes[genes %in% rownames(expr_matrix)]
  })
  
  # 计算每个通路的平均表达值
  pathway_activity <- lapply(pathway_list, function(genes) {
    colMeans(expr_matrix[genes, , drop = FALSE])
  })
  
  # 转换为数据框
  pathway_activity_df <- as.data.frame(do.call(rbind, pathway_activity))
  rownames(pathway_activity_df) <- names(pathway_list)
  
  # 为每个通路计算模块得分
  seurat <- AddModuleScore(
    object = seurat,
    features = pathway_list,
    name = c("Pathway1", "Pathway2")
  )
  
  # 查看 meta.data 中的新列名
  colnames(seurat@meta.data)
  
  # 假设输出的列名为 "Pathway11" 和 "Pathway21"
  head(seurat@meta.data[, c("Pathway11", "Pathway22")])
  
  
  library(ggplot2)
  
  # 添加分组信息
  meta_data <- seurat@meta.data
  
  # 绘制箱线图
  p <- ggplot(meta_data, aes(x = condition, y = Pathway11, fill = condition)) +
    geom_violin(trim = FALSE, alpha = 0.3, scale = "width", linewidth = 0.35) + # 增加小提琴图，透明度调低
    geom_boxplot(aes(fill = condition),
                 alpha = 0.5,
                 linewidth = 0.5,
                 notch = FALSE,
                 width = 0.5,
                 outlier.size = 0.05,
                 outlier.color = "black"
    ) +
    geom_jitter(mapping = aes(colour = condition),
                shape = 16,
                position = position_jitter(0.2),
                size = 0.3,
                alpha = 0.8
    ) +
    geom_signif(
      comparisons = list(c("Healthy", "Disease")), 
      map_signif_level = T, 
      test = "t.test", # 或 "wilcox.test" 视具体情况
      step_increase = 1,
      tip_length = 0,
      textsize = 3.5,
      size = 0.35
    ) +
    stat_summary(fun = mean, geom = "point", color = "blue", size = 2) +
    stat_summary(fun = median, geom = "point", color = "red", size = 2) +
    labs(
      title = "Cellular senescence Activity by Condition", 
      x = "Condition", 
      y = "Cellular senescence Activity Score"
    ) +
    theme(
      # panel.grid.major = element_line(color = "grey90", size = 0.35),
      # panel.border = element_rect(color = "black", fill = NA, size = 0.35),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line(linewidth = 0.35),
      axis.ticks = element_line(linewidth = 0.35, colour = "grey0"),
      text = element_text(size = 8, colour = "grey0"),
      axis.text.x = element_text(size = 8,  colour = "grey0"),
      axis.text.y = element_text(size = 8, colour = "grey0"),
      axis.title.y = element_text(size = 8, colour = "grey0"),
      axis.title.x = element_blank(),
      plot.title = element_text(size = 8, fAEC = "bold"),
      legend.position = "none"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
    # scale_y_log10(expand = expansion(mult = c(0, 0.1))) +
    scale_fill_manual(values = c("Healthy" = "#87cEEB", "Disease" = "#FF5722")) +
    scale_colour_manual(values = c("Healthy" = "#87cEEB", "Disease" = "#FF5722"))
  # 显示图
  print(p)
  # 保存图片,多个柱子，宽度加0.5
  ggsave(paste0("Cellular senescence Activity Score", ".tiff"), plot = p, height = 3, width = 3,dpi = 300)
  ggsave(paste0("Cellular senescence Activity Score", ".pdf"), plot = p, height = 3, width = 3,dpi = 300)
  
  
  
  
  t.test(meta_data$Pathway11 ~ meta_data$condition)  # Pathway1 活性
  t.test(meta_data$Pathway22 ~ meta_data$condition)  # Pathway2 活性
  cor(meta_data$Pathway11, meta_data$Pathway22, method = "pearson")
  
  
  
  ggplot(meta_data, aes(x = Pathway11, y = Pathway22)) +
    geom_point(alpha = 0.6) +
    stat_smooth(method = "lm", color = "blue") +
    labs(title = "Correlation between Pathway1 and Pathway2",
         x = "Pathway1 Activity", y = "Pathway2 Activity") +
    stat_cor(method = "pearson")
  
  
  
  
  
  # 提取 AT2 细胞的表达矩阵----
  at2_expr_matrix <- GetAssayData(seurat, layer = "data")
  
  # 设定阈值筛选高表达基因
  at2_high_genes <- rownames(at2_expr_matrix)[rowMeans(at2_expr_matrix) > 0.2]  # 示例阈值
  
  # 筛选 AT2 高表达基因与通路基因的交集
  pathway1_genes_filtered <- intersect(pathway1_genes, at2_high_genes)
  pathway2_genes_filtered <- intersect(pathway2_genes, at2_high_genes)
  
  # 分组标注（假设有 Healthy 和 Disease）
  Idents(seurat) <- "condition"
  
  # 计算分组间的差异表达
  de_genes <- FindMarkers(seurat, ident.1 = "Disease", ident.2 = "Healthy")
  
  # 筛选显著差异基因（假设 p < 0.05）
  de_genes_filtered <- rownames(de_genes[de_genes$p_val_adj < 0.05, ])
  
  # 在这些基因中提取 GSEA 显著通路
  gsea_result <- gseGO(geneList = de_genes_filtered, OrgDb = org.Hs.eg.db, keyType = "SYMBOL")
  
  # 热图绘制
  p3 <- heatplot(tmp, foldChange=geneList, showCategory=5,label_format = 30)
  p3
  # 查看通路
  selected_pathway <- geneSetIDs
  for (pathway.id in selected_pathway ){
    print(pathway.id)
    pathview(gene.data = geneList,
             pathway.id = pathway.id,
             species = "mmu",
             limit         = list(gene=max(abs(geneList)), cpd=1))
    pathview(gene.data     = geneList,
             pathway.id    = pathway.id,
             species       = 'mmu' ,      # 人类hsa 小鼠mmu 
             kegg.native   = F,# TRUE为png，F为pdf
             new.signature = T, #pdf是否显示pathway标注
             limit         = list(gene=max(abs(geneList)), cpd=1)#图例color bar范围调整 
    )
  }
  
  # 上下调分开----
  dotplot(
    gesKK,
    font.size = 8,
    showCategory = 10,
    split = ".sign")+fAECt_grid(.~.sign)#上下调通路分开
  
  cnetplot(gesKK1,foldchange = geneList,
           circular = TRUE,
           color.params = list(edge = TRUE,
                               category = "#E5C494",
                               gene = "#B3B3B3"),
           showCategory = 5,
           cex.params = list(gene_label = 0.8,
                             category_label = 0.8,
                             category_node = 0.8,
                             gene_node = 0.8)
  )
  ggsave("KEGG_res_circle.pdf",width=13,height=9)
  dev.off()
  # GSEA_kk, Dotplot 自定义气泡图----
  library(tidyverse)
  # desc <- sub("-.*", "", gseKK_res$Description)  # 提取连字符前的部分
  Dotplot <- function(gseKK_res, title="") {
    dotplot <- ggplot(cbind(gseKK_res, Order = nrow(gseKK_res):1)) +
      geom_point(mapping = aes(x = -log10(p.adjust), y = Order, 
                               size = setSize, fill = enrichmentScore),
                 shape = 21) + 
      scale_fill_gradientn(colours = c("blue", "red")) + #自定义配色
      scale_y_continuous(position = "left", 
                         breaks = 1:nrow(gseKK_res), 
                         labels = Hmisc::capitalize(rev( gseKK_res$Description))) +
      # scale_x_continuous(breaks = c(3, 4,5,6),
      #                   #breaks = seq(0, xmax+5, 5),
      #                   limits = c(3,6),
      #                   expand = expansion(mult = c(.05, .05))) + #两边留空
      labs(x = "-Log10(p.adjust)", y = NULL) +
      guides(size = guide_legend(title = "setSize"),
             fill = guide_colorbar(title = "enrichmentScore")) +
      ggtitle(title) +
      theme_bw() +
      theme(panel.grid =element_blank(),
            panel.grid.major = element_line(color = "gray", linetype = "dashed"),
            panel.border = element_rect(color = "black", linewidth = 1),
            axis.text = element_text(size = rel(1)),
            title = element_text(size = rel(1))) #去除网格线
    # dotplot %>% ggplotGrob()#  %>% cowplot::plot_grid()
    # # fname=paste0(resultdir, '/', filemark,'.pdf')
    # # ggsave(fname, width = 8, height = 10)
    return(dotplot)
  }
  range(round(-log10(gseKK_res$p.adjust)))
  p_df <- gseKK_res[1:20,]
  Dotplot(p_df, title="GSEA-KEGG Enrich")
  ggsave(file = "Dotplot_GSEA-kegg_top20_enrichmentscore.pdf", width = 11, height = 7)
  #gh <- ggplotGrob(dotplot)
  #gd <- ggplotGrob(dotplotk)
  #cowplot::plot_grid(gh, gd, rel_widths = c(1.2, 1))
  
  sjPlot::tab_gseGO(title = '',gseKK_res_order[,2:8]) %>% #使用三线表查看结果

}

#展示多条通路
library(ggplotify)
library(msigdbr)
# GSEA,msigdbr，GO富集和clusterProfiler直接分析的GSEA_GO一致----
genesets_M5 <- clusterProfiler::read.gmt("m8.all.v2023.2.Mm.entrez.gmt")
msigdbr_species() #没更新
genesets <- msigdbr(species = "Homo sapiens",
                    category = "H") %>%
  dplyr::select(gs_name, entrez_gene)
gsea_res <- GSEA(geneList, 
                 TERM2GENE = genesets,
                 minGSSize = 1,
                 maxGSSize = 500,
                 pvalueCutoff = 0.9,
                 pAdjustMethod = "BH")
# gsea_res[[gsea_res$ID[[1]]]]
tmp <- setReadable(gsea_res, "org.Hs.eg.db", "ENTREZID")
gsea_res <- tmp@result %>% dplyr::arrange(desc(enrichmentScore))
head(gsea_res$Description,3)
writexl::write_xlsx(gsea_res,"msigdbr_GO_res.xlsx", format_headers = T,col_names = T)

writexl::write_xlsx(gsea_res, "msigdbr_KEGG_res.xlsx", format_headers = T,col_names = T)
# 绘制GSEA 图
# 排序
gsea_res_order <- gsea_res %>% dplyr::arrange(desc(enrichmentScore))
# 绘图代码
geneSetIDs <- c("HALLMARK_P53_PATHWAY")
if(T){
  p <- gseaplot2(tmp, geneSetID = geneSetIDs,
                 title = "",
                 base_size = 8,
                 pvalue_table = T,
                 rel_heights = c(1.5,0.5,1),
                 color = colorspAEC::rainbow_hcl(3),
                 
                 ES_geom = "line")
  p[[1]] <- p[[1]] +  
    # scale_color_viridis_d(labels = tmp$Description)+ #颜色参数
    geom_hline(yintercept = 0,
               color = "grey75", 
               linewidth = 0.35, 
               linetype = 2)+
    theme(legend.position = "top", 
          legend.direction = "horizontal")
  # p[[2]] <- p[[2]] + 
  # scale_color_viridis_d(labels = tmp$Description)#颜色参数
  p[[3]] <- p[[3]] +
    geom_hline(yintercept = 0, color = "steelblue", 
               linewidth = 0.35,
               linetype = 2)
  p
  ggsave("GSEA_P53.pdf",width=8,height=6)
  ggsave("GSEA_P53.tiff",width=8,height=6)
}

# 热图绘制
p3 <- heatplot(tmp, foldChange = geneList, showCategory = 6,label_format = 30)
#GSEA结果保存####
write.csv(gsea_res_order, file = "GSEA_results.csv", row.names = T)####
sjPlot::tab_df(title = '',gsea_res_order[1:10,2:8])#使用三线表查看结果


##查看与选择所需通路####
gsekk_res <- gesKK@result
i = 1
select_pathway <- gesKK@result$ID[1] #选择所需通路的ID号
pathview(gene.data     = geneList,
         pathway.id    = "mmu04657",
         out.suffix = "map",
         species       = 'mmu' ,      # 人类hsa 小鼠mmu 
         kegg.native   = T,# TRUE输出完整pathway的png文件，F输出基因列表的pdf文件 
         new.signature = F, #pdf是否显示pathway标注
         limit         = list(gene=2.5, cpd=1)
)
# 然后循环绘图####
###设置工作目录到想要的地方
dir.create("C:/WJ/Experiment Data-WJ/RNASeq/pathwayview_out")
setwd("C:/WJ/Experiment Data-WJ/RNASeq/pathwayview_out")

for (pathway.id in gesKK@result$ID[1:5]){
  print(pathway.id)
  pathview(gene.data = geneList,
           pathway.id = pathway.id,
           species = "mmu",
           limit         = list(gene=max(abs(geneList)), cpd=1))
  pathview(gene.data     = geneList,
           pathway.id    = pathway.id,
           species       = 'mmu' ,      # 人类hsa 小鼠mmu 
           kegg.native   = F,# TRUE为png，F为pdf
           new.signature = T, #pdf是否显示pathway标注
           limit         = list(gene=max(abs(geneList)), cpd=1)#图例color bar范围调整 
  )
}


setwd("C:/WJ/Experiment Data-WJ/RNASeq")

## 基因概念网####
pdf(file="基因概念网络gsekk.pdf",width=13,height=8)
geom_text_repel(max.overlaps = 50000, nudge_x = 1, nudge_y = 1)
cnetplot(gesKK1,foldchange = geneList,
         circular = TRUE, 
         color.params = list(edge = TRUE),
         showCategory = 5)
dev.off()
heatplot(gesKK1, foldChange=geneList)

