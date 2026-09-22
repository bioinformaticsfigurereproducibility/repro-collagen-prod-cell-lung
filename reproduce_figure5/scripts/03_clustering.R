# ==============================================================================
# Integration and Clustering of _All + _Lin Dataset
#   "Performed graph-based clustering of the MNN reduced data using FindClusters
#   function of Seurat (resolution= 0.3, dims.use= 1:19)"
# ==============================================================================

# ---- 0. Load Libraries & Checkpoint ------------------------------------------
library(SeuratWrappers)
library(Seurat)
library(glmGamPoi)
library(ggplot2)

merged_all_sce <- readRDS("reproduce_figure5/data/filtered_sce_scran_qc.rds")

# ---- 1. Convert into Seurat Object -------------------------------------------
merged_all <- as.Seurat(merged_all_sce, counts = "counts", data = NULL)

# ---- 2. Ensure Seurat Object has correct metadata ----------------------------
merged_all$orig.ident <- merged_all$sample
Idents(merged_all) <- "orig.ident"
table(merged_all$condition)
merged_all$percent.mt <- merged_all$subsets_Mito_percent

Assays(merged_all)
merged_all <- RenameAssays(merged_all, originalexp = "RNA")

# ---- 3. Normalize, Scale and Find Variable Features through SCTransform() ----
# To reduce computational capacity, run SCTransform per-sample split layers
merged_all[["RNA"]] <- split(merged_all[["RNA"]], f = merged_all$orig.ident)
merged_all <- SCTransform(merged_all, vars.to.regress = "percent.mt", vst.flavor = "v2")

merged_all <- RunPCA(merged_all)

all_elbow_plot <- ElbowPlot(merged_all, ndims = 50)
ggsave(
  filename = "reproduce_figure5/figures/elbow_plot.png",
  plot = elbow_plot,
  width = 8,
  height = 6
)

# ---- 4. Checkpoint -----------------------------------------------------------
saveRDS(merged_all, file = "reproduce_figure5/data/merged_all_sct.rds")

# ---- 5. MNN Integration ------------------------------------------------------
# MNN Integration requires multiple batches
merged_all[["SCT"]] <- split(merged_all[["SCT"]], f = merged_all$orig.ident)
# Changes from SCTAssay to Assay5
Layers(merged_all[["SCT"]])

merged_all <- IntegrateLayers(
  object = merged_all,
  method = FastMNNIntegration,
  orig.reduction = "pca",
  new.reduction = "mnn",
  verbose = FALSE
)

# Check for completion
Reductions(merged_all)

# ---- 6. Standardized Workflow ------------------------------------------------
merged_all <- FindNeighbors(merged_all, reduction = "mnn", dims = 1:19)
merged_all <- FindClusters(merged_all, resolution = 0.3, graph.name = "RNA_snn")
merged_all <- RunUMAP(merged_all, reduction = "mnn", dims = 1:19)

# ---- 7. Visualization --------------------------------------------------------
all_dim_plot <- DimPlot(merged_all, reduction = "umap", label = TRUE)
ggsave(
  filename = "reproduce_figure5/figures/all_dim_plot.png",
  plot = dim_plot,
  width = 8,
  height = 6
)
all_dim_plot_grouped_condition <- DimPlot(merged_all, reduction = "umap", group.by = "condition")
ggsave(
  filename = "reproduce_figure5/figures/all_dim_plot_grouped_condition.png",
  plot = dim_plot_grouped_condition,
  width = 8,
  height = 6
)
all_dim_plot_split_condition <- DimPlot(merged_all, reduction = "umap", split.by = "condition")
ggsave(
  filename = "reproduce_figure5/figures/all_dim_plot_split_condition.png",
  plot = dim_plot_split_condition,
  width = 8,
  height = 6
)

# ---- 8. Checkpoint -----------------------------------------------------------
saveRDS(merged_all, file = "reproduce_figure5/data/merged_all_integrated.rds")
