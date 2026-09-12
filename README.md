# scRNA-seq Analysis: Breast Cancer Circulating Tumour Cells (CTCs)

## Overview
This project reproduces the key findings of Jordan et al. (2016) using 
single cell RNA sequencing data from circulating breast cancer cells.
The goal was to identify a HER2+ subpopulation hiding within clinically 
classified HER2- breast cancer patients.

## Key Finding
We successfully identified two distinct cell populations:
- **HER2+ cells** (38 cells) — expressing high ERBB2, epithelial markers
- **HER2- cells** (23 cells) — expressing immune related markers

## Dataset
- **Source:** ArrayExpress E-GEOD-75367
- **Link:** https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-GEOD-75367
- **Cells:** 74 circulating tumour cells from breast cancer patients
- **Paper:** Jordan et al. Nature 537, 102-106 (2016)

## Methods
| Step | Method |
|---|---|
| Quality Control | Filtered cells by gene count and MT% |
| Normalisation | LogNormalize (scale factor 10,000) |
| Feature Selection | Top 3,000 variable genes (VST) |
| Dimensionality Reduction | PCA (15 PCs) |
| Clustering | Louvain algorithm (resolution 0.8) |
| Visualisation | UMAP and tSNE |
| Marker Discovery | FindAllMarkers (Wilcoxon test) |

## Tools
- R version 4.5.1
- Seurat v5
- clusterProfiler
- org.Hs.eg.db

## Key Marker Genes Found

### HER2+ Markers
| Gene | Function | Fold Change |
|---|---|---|
| KRT19 | Epithelial cancer marker | 286x |
| AZGP1 | Breast cancer prognosis | 194x |
| APOD | Lipid transport | 161x |
| CRABP2 | Cell differentiation | 160x |
| KRT7 | Epithelial marker | 154x |

### HER2- Markers
| Gene | Function | Fold Change |
|---|---|---|
| PLEK | Immune signalling | 136x |
| STX11 | Cell secretion | 99x |
| TM6SF1 | Fat metabolism | 50x |
| SUSD3 | Cell surface protein | 49x |
| LGALS12 | Cell death signals | 48x |

## Output Plots
| Plot | Description |
|---|---|
| 01_QC_before_filtering | Cell quality before filtering |
| 02_QC_after_filtering | Cell quality after filtering |
| 03_variable_features | Most variable genes |
| 04_PCA_loadings | Genes driving each PC |
| 05_PCA_plot | Cells in PCA space |
| 06_elbow_plot | PC selection guide |
| 07_UMAP_clusters | Clusters on UMAP |
| 08_UMAP_ERBB2 | ERBB2 expression on UMAP |
| 09_ERBB2_violin | ERBB2 per cluster |
| 10_UMAP_final_labelled | Final HER2+/HER2- map |
| 11_marker_heatmap | Marker gene heatmap |
| 12_dotplot_markers | Summary dotplot |
| 13_PCA_labelled | PCA with labels |
| 14_PCA_ERBB2 | ERBB2 on PCA |

## Citation
Jordan, Nicole Vincent et al. 
"HER2 expression identifies dynamic functional states within 
circulating breast cancer cells." 
Nature vol. 537,7618 (2016): 102-106.

## Author
Akingbade Boluwatife Samuel
GitHub: https://github.com/Akingtom
