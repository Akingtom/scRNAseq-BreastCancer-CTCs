# ==============================================================================
# Single-Cell Transcriptomic Analysis of Breast Cancer CTCs (E-GEOD-75367)
# Pipeline: QC -> Normalization -> DimRedux -> Differential Expression ->
#           Functional Enrichment (GO/KEGG) -> Clinical Metadata Alignment
# ==============================================================================

# 1. Load Required Libraries ---------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

# Set working directory/path
base_path <- "."

# 2. Ingest Data & Construct Count Matrix ---------------------------------------
counts_raw <- readMM(file.path(base_path, "E-GEOD-75367-quantification-raw-files.counts.mtx"))
barcodes   <- readLines(file.path(base_path, "E-GEOD-75367-quantification-raw-files.counts.mtx_cols"))
features   <- readLines(file.path(base_path, "E-GEOD-75367-quantification-raw-files.counts.mtx_rows"))

colnames(counts_raw) <- barcodes
rownames(counts_raw) <- features

# Map Ensembl Gene IDs to HGNC Symbols
gene_map <- mapIds(
  org.Hs.eg.db,
  keys = features,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# Resolve unmapped and duplicated gene symbols
clean_gene_names <- ifelse(is.na(gene_map), features, gene_map)
clean_gene_names <- make.unique(clean_gene_names)
rownames(counts_raw) <- clean_gene_names

# Create Seurat Object
BCdata <- CreateSeuratObject(counts = counts_raw, project = "BC_CTCs", min.cells = 3)

# 3. Quality Control Filtering -------------------------------------------------
BCdata[["percent.mt"]] <- PercentageFeatureSet(BCdata, pattern = "^MT-")

# Filter low-quality cells / artifacts
BCdata <- subset(
  BCdata,
  subset = nFeature_RNA >= 200 & 
    nCount_RNA < 5000000 & 
    percent.mt < 20
)

# 4. Normalization, Feature Selection & Dimensionality Reduction ---------------
BCdata <- NormalizeData(BCdata, normalization.method = "LogNormalize", scale.factor = 10000)
BCdata <- FindVariableFeatures(BCdata, selection.method = "vst", nfeatures = 3000)
BCdata <- ScaleData(BCdata, features = rownames(BCdata))
BCdata <- RunPCA(BCdata, features = VariableFeatures(object = BCdata), npcs = 15)

# Graph-based Clustering & UMAP
BCdata <- FindNeighbors(BCdata, dims = 1:15)
BCdata <- FindClusters(BCdata, resolution = 0.8)
BCdata <- RunUMAP(BCdata, dims = 1:15)

# 5. Cell Phenotyping & Subtype Annotation ------------------------------------
# Identify phenotypic markers across clusters
markers <- FindAllMarkers(
  BCdata,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.5
)

# Annotate clusters based on ERBB2 / epithelial vs platelet cloaking signatures
# (Cluster 0: Epithelial/HER2+, Cluster 1: Platelet-associated HER2-)
new_cluster_ids <- c("HER2+", "HER2-")
names(new_cluster_ids) <- levels(BCdata)
BCdata <- RenameIdents(BCdata, new_cluster_ids)
BCdata$phenotype <- Idents(BCdata)

# 6. Functional Enrichment Analysis (GO & KEGG) --------------------------------
# Convert marker symbols to Entrez IDs for over-representation tests
cluster0_genes <- markers %>% filter(cluster == levels(markers$cluster)[1]) %>% pull(gene)
cluster1_genes <- markers %>% filter(cluster == levels(markers$cluster)[2]) %>% pull(gene)

entrez_c0 <- mapIds(org.Hs.eg.db, keys = cluster0_genes, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
entrez_c1 <- mapIds(org.Hs.eg.db, keys = cluster1_genes, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")

entrez_c0 <- na.omit(entrez_c0)
entrez_c1 <- na.omit(entrez_c1)

# GO Biological Process
go_her2_pos <- enrichGO(gene = entrez_c0, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05)
go_her2_neg <- enrichGO(gene = entrez_c1, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05)

# KEGG Pathways
kegg_her2_pos <- enrichKEGG(gene = entrez_c0, organism = "hsa", pvalueCutoff = 0.05)
kegg_her2_neg <- enrichKEGG(gene = entrez_c1, organism = "hsa", pvalueCutoff = 0.05)

# Generate and save pathway visualization
p_kegg_pos <- dotplot(kegg_her2_pos, showCategory = 8) + 
  scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +
  ggtitle("KEGG (HER2+ Markers)")

p_kegg_neg <- dotplot(kegg_her2_neg, showCategory = 8) + 
  scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +
  ggtitle("KEGG (HER2- Markers)")

ggsave("kegg_enrichment_dotplot.png", plot = p_kegg_pos / p_kegg_neg, width = 9, height = 7)

# 7. Clinical SDRF Metadata Mapping & Batch Assessment ------------------------
sdrf <- read_delim(file.path(base_path, "E-GEOD-75367.sdrf.txt"), delim = "\t", trim_ws = TRUE)

# Strip punctuation/spaces for exact matching
clean_sdrf_keys <- gsub("[^[:alnum:]]", "", tolower(sdrf$`Comment [Sample_title]`))
clean_cell_keys <- gsub("[^[:alnum:]]", "", tolower(colnames(BCdata)))

matched_idx  <- match(clean_cell_keys, clean_sdrf_keys)
sdrf_matched <- sdrf[matched_idx, ]

# Map individual identifiers with barcode fallback
BCdata$patient_id <- sdrf_matched$`Characteristics [individual]`
fallback_ids <- sub("^(BRx[-_][0-9]+).*", "\\1", colnames(BCdata))
fallback_ids <- gsub("-", "_", fallback_ids)
BCdata$patient_id[is.na(BCdata$patient_id)] <- fallback_ids[is.na(BCdata$patient_id)]

# Final UMAP assessment: Patient distribution vs Biological clusters
p_patient <- DimPlot(BCdata, reduction = "umap", group.by = "patient_id") + ggtitle("Cells by Patient ID")
p_cluster <- DimPlot(BCdata, reduction = "umap") + ggtitle("Transcriptomic Clusters")

ggsave("umap_patient_vs_cluster.png", plot = p_patient + p_cluster, width = 11, height = 5)

# 8. Export Results ------------------------------------------------------------
write.csv(markers, "HER2_cluster_markers.csv", row.names = FALSE)
saveRDS(BCdata, "BCdata_final_with_patients.rds")

message("Workflow execution complete. Results and figures saved.")