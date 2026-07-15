# GSE193816 — Single-cell RNA-seq Analysis

This repository contains the analysis pipeline for the **GSE193816** scRNA-seq dataset, profiling immune and epithelial cells in asthma patients and healthy controls.

## Dataset

GSE193816: Single-cell transcriptomic profiling of induced sputum from asthma patients and healthy controls.

- **Samples**: Ag (antigen-stimulated), Pre (pre-vaccination), Dil (diluent — excluded)
- **Conditions**: Asthma vs Control
- **Input format**: h5ad (converted to Seurat via `sceasy`)

## Pipeline Structure

```
pipeline/
├── 01_load_and_qc.R          # Data loading, QC filtering, doublet removal
├── 02_cluster_and_markers.R  # Clustering, UMAP, marker gene discovery
├── 03_cell_annotation.R      # Manual cell type annotation + visualization
├── 04_differential_expression.R  # DE analysis, violin/box plots
├── 05_enrichment_analysis.R  # GSEA / GO / KEGG pathway enrichment
├── R/
│   └── utils.R               # Shared plotting and annotation functions
├── results/                  # Output directory (plots, CSVs, .qs files)
├── .gitignore
└── README.md
```

## Requirements

### Core packages
```r
install.packages(c("Seurat", "tidyverse", "qs", "RColorBrewer", "patchwork",
                   "ggsignif", "BiocManager"))
```

### Bioconductor packages
```r
BiocManager::install(c("sceasy", "scDblFinder", "SingleR", "celldex",
                        "clusterProfiler", "enrichplot", "org.Hs.eg.db",
                        "presto", "BiocParallel"))
```

### External tool
The pipeline references `seurat_toolbox.R` for QC functions (`basic_qc`, `filter_by_gene_umi`, `subType`). Place this file in `pipeline/R/` or update the `source()` path in `01_load_and_qc.R`.

## Usage

1. Place `GSE193816_all_data_raw_counts.h5ad` in `pipeline/data/`
2. Run scripts in order:

```bash
Rscript pipeline/01_load_and_qc.R
Rscript pipeline/02_cluster_and_markers.R
Rscript pipeline/03_cell_annotation.R
Rscript pipeline/04_differential_expression.R
Rscript pipeline/05_enrichment_analysis.R
```

Or source them interactively from RStudio.

## Output

- `results/seurat_anno.qs` — Annotated Seurat object
- `results/global_markers.csv` — Marker genes per cluster
- `results/umap_annotated.png` — UMAP with cell type labels
- `results/dotplot_markers.png` — Marker gene dotplot
- `results/GSEA_*` — Pathway enrichment results

## Notes

- The `scripts/01_load_and_qc.R` requires `seurat_toolbox.R` (external file with custom QC functions).
- Cluster-to-cell-type mappings in `03_cell_annotation.R` should be validated and adjusted for your specific clustering results.
