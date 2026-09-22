# ==============================================================================
# Reproduce Figure 5c-g: Subset of COL1A1 genes within the _Lin dataset
#   "Visualization of the clusters on a 2D map was performed with UMAP 
#   (RunUMAP function of Seurat, dims.use= 1:19)"
#
#   Fig 5c: UMAP by condition
#   Fig 5d: UMAP by cluster
#   Fig 5e/5g: FeaturePlots of Marker genes
#   Fig 5f: COL1A1/CTHRC1/ACTA2 split by condition
# ==============================================================================

# ---- 0. Load Libraries and Checkpoint ----------------------------------------
library(SeuratWrappers)
library(Seurat)
library(glmGamPoi)
library(ggplot2)
library(patchwork)

merged_all <- readRDS(file = "reproduce_figure5/data/merged_all_final_annotated.rds")

# ---- 1. Filter to only Lin-only Samples --------------------------------------
lin_samples <- grep("_Lin$", unique(merged_all$orig.ident), value = TRUE)
merged_lin <- subset(merged_all, subset = orig.ident %in% lin_samples)

# ---- 2. Join Layer to Subset Dataset ----------------------------------------
merged_lin <- JoinLayers(merged_lin, assay = "RNA")
DefaultAssay(merged_lin) <- "RNA"

# ---- 3. Subset to COL1A1 genes -----------------------------------------------
merged_lin_col1a1 <- subset(merged_lin, subset = COL1A1 > 0)

# Paper stated 48,587 cells, after subsetting only 25,782 cells
# Could be due to pipeline differences (all NML samples showed lower COL1A1 levels, consistent with biology)
# COL1A1 encodes type 1 collagen which is found in fibrotic scar tissue, producing scar-forming collagen
ncol(merged_lin_col1a1)
table(merged_lin_col1a1$orig.ident)

# ---- 4. Re-conduct SCTranform focused on COL1A1 population -------------------
merged_lin_col1a1[["RNA"]] <- split(merged_lin_col1a1[["RNA"]], f = merged_lin_col1a1$orig.ident)
merged_lin_col1a1 <- SCTransform(merged_lin_col1a1, vars.to.regress = "percent.mt", vst.flavor = "v2")

# --- 5. Pre-processing --------------------------------------------------------
merged_lin_col1a1 <- RunPCA(merged_lin_col1a1)
col1a1_elbow_plot <- ElbowPlot(merged_lin_col1a1, ndims = 50)
ggsave(
  filename = "reproduce_figure5/figures/col1a1_elbow_plot.png",
  plot = col1a1_elbow_plot,
  height = 6,
  width = 8
)

# ---- 6. Checkpoint -----------------------------------------------------------
saveRDS(merged_lin_col1a1, "reproduce_figure5/data/merged_lin_col1a1_sct.rds")

# ---- 7. MNN Integration ------------------------------------------------------
merged_lin_col1a1[["SCT"]] <- split(merged_lin_col1a1[["SCT"]], f = merged_lin_col1a1$orig.ident)

merged_lin_col1a1 <- IntegrateLayers(
  object = merged_lin_col1a1,
  method = FastMNNIntegration,
  orig.reduction = "pca",
  new.reduction = "mnn",
  verbose = FALSE
)

# ---- 8. Standardized Workflow ------------------------------------------------
merged_lin_col1a1 <- FindNeighbors(merged_lin_col1a1, reduction = "mnn", dims = 1:19)
merged_lin_col1a1 <- FindClusters(merged_lin_col1a1, resolution = 0.3, graph.name = "RNA_snn")
merged_lin_col1a1 <- RunUMAP(merged_lin_col1a1, reduction = "mnn", dims = 1:19)

# ---- 9. Visualization --------------------------------------------------------
## --- 9a. Figure 5c -----------------------------------------------------------
col1a1_dimplot_grouped_condition <- DimPlot(merged_lin_col1a1, 
                                            reduction = "umap", 
                                            group.by = "condition") + ggtitle("COL1A1+ Cells")
ggsave(
  filename = "reproduce_figure5/figures/fig5c.png",
  plot = col1a1_dimplot_grouped_condition,
  height = 6,
  width = 8
)
## --- 9b. Figure 5d -----------------------------------------------------------
col1a1_dimplot <- DimPlot(merged_lin_col1a1, reduction = "umap", label = TRUE)
ggsave(
  filename = "reproduce_figure5/figures/fig5d.png",
  plot = col1a1_dimplot,
  height = 6,
  width = 8
)

## ---- 9c. Figure 5e & Figure 5f ----------------------------------------------------------
genes_e <- c("PDGFRA", "MYH11", "RGS5", "NPNT", "PI16", "FGF18",
             "PDGFRB", "MYL9", "MCAM", "CES1", "DCN", "WIF1")
genes_g <- c("TNC", "POSTN", "COL3A1")

# To have same scaled expression through each panel
all_genes <- unique(c(genes_e, genes_g))
global_max <- max(FetchData(merged_lin_col1a1, vars = all_genes))
fig5e <- FeaturePlot(merged_lin_col1a1, features = genes_e, reduction = "umap",
                     order = TRUE, ncol = 6) & NoAxes() &
  scale_color_gradientn(colors = c("lightgrey", "blue"), limits = c(0, global_max))

fig5g <- FeaturePlot(merged_lin_col1a1, features = genes_g, reduction = "umap",
                     order = TRUE, ncol = 3) & NoAxes() &
  scale_color_gradientn(colors = c("lightgrey", "blue"), limits = c(0, global_max))

ggsave(
  filename = "reproduce_figure5/figures/fig5e.png", 
  plot = fig5e, 
  width = 18, 
  height = 6)

ggsave(
  filename = "reproduce_figure5/figures/fig5g.png", 
  plot = fig5g, 
  width = 9, 
  height = 3)

## ---- 9e. Figure 5f ----------------------------------------------------------
# Same order of condition as the paper (NML, IPF, SCD)
merged_lin_col1a1$condition <- factor(merged_lin_col1a1$condition, levels = c("NML", "IPF", "SCD"))
fig5f <- FeaturePlot(
  merged_lin_col1a1,
  features   = c("COL1A1", "CTHRC1", "ACTA2"),
  split.by   = "condition",
  reduction  = "umap",
  order      = TRUE,
  cols       = c("lightgrey", "red"),
  keep.scale = "feature"
) & NoAxes()

ggsave(
  filename = "reproduce_figure5/figures/fig5f.png", 
  plot = fig5f, 
  width = 12, 
  height = 9)

# ---- 10. In-Depth Analyzing Cluster 7 & 8 ------------------------------------
merged_lin_col1a1_2 <- JoinLayers(merged_lin_col1a1, assay = "SCT")

cluster8_markers <- FindMarkers(merged_lin_col1a1_2, ident.1 = "8", ident.2 = c("0", "5"))
head(cluster8_markers[order(-cluster8_markers$avg_log2FC), ], 20)
FeaturePlot(merged_lin_col1a1_2, features = c("WT1", "MSLN", "UPK3B"), reduction = "umap")
DotPlot(merged_lin_col1a1_2, features = c("WT1", "MSLN", "UPK3B"), group.by = "seurat_clusters") + RotatedAxis()
saveRDS(merged_lin_col1a1, "reproduce_figure5/data/merged_lin_col1a1_clustered.rds")