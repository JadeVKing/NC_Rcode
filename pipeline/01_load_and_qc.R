# =============================================================================
# GSE193816 scRNA-seq: Data loading and QC
# =============================================================================
# Requirements:
#   - Input: GSE193816_all_data_raw_counts.h5ad  (in data/)
#   - seurat_toolbox.R  (provide path via source() or place in R/)
#   - Packages: Seurat, sceasy, scDblFinder, harmony, qs
# =============================================================================

library(Seurat)
library(sceasy)
library(scDblFinder)
library(BiocParallel)
library(qs)

# ---- Configuration ----
data_dir <- "data"
results_dir <- "results"
input_h5ad <- file.path(data_dir, "GSE193816_all_data_raw_counts.h5ad")

# Source QC toolbox (adjust path as needed)
source("R/seurat_toolbox.R")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 1. Convert h5ad to Seurat ----
register(SnowParam(workers = 4, progressbar = TRUE))

sceasy::convertFormat(
  obj = input_h5ad,
  from = "anndata",
  to = "seurat",
  outFile = file.path(data_dir, "scRNA.rds")
)

seurat <- readRDS(file.path(data_dir, "scRNA.rds"))
DefaultAssay(seurat) <- "RNA"

# ---- 2. Inspect metadata ----
colnames(seurat@meta.data)
head(seurat@meta.data)
table(seurat$sample)

# ---- 3. Subset to Ag and Pre samples (exclude Dil) ----
seurat <- subset(seurat, subset = sample %in% c("Ag", "Pre"))
seurat$sample <- droplevels(seurat$sample)
table(seurat$sample)

# ---- 4. QC ----
seurat <- basic_qc(seurat)
seurat <- filter_by_gene_umi(seurat)
seurat <- subType(seurat)

# ---- 5. Doublet detection ----
table(seurat$id)
sce <- as.SingleCellExperiment(seurat)
sce <- scDblFinder(sce, samples = seurat$Channel)
seurat$scDblFinder.class <- sce$scDblFinder.class
rm(sce); gc()

# Visualize doublet calls
CellDimPlot(seurat, group.by = "scDblFinder.class", reduction = "umap")
table(seurat$scDblFinder.class)

# ---- 6. Remove doublets and re-cluster ----
seurat <- subset(seurat, subset = scDblFinder.class == "singlet")
DefaultAssay(seurat) <- "RNA"
seurat <- subType(seurat,
  res = 0.8,
  sample_cells = 200000,
  use_harmony = FALSE,
  harmony_group = "sample_id"
)

seurat <- FindNeighbors(seurat, dims = 1:15)
seurat <- FindClusters(seurat, resolution = 0.8)
seurat <- RunUMAP(seurat, dims = 1:16, n.neighbors = 30)

# ---- 7. Save clean object ----
qsave(seurat, file.path(results_dir, "seurat_clean_Ag&Pre.qs"))

message("QC and filtering complete. Output: ", file.path(results_dir, "seurat_clean_Ag&Pre.qs"))
