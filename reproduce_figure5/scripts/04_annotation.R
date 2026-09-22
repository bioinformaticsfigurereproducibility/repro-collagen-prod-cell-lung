# ==============================================================================
# Annotation - Automatic and Manual Annotations
#   Usage of Azimuth (Automatics) vs Manual (CellMarker 3.0)
#   Compared side by side  and disagreements were resolved by checking marker
#   genes directly. Final Annotation reflects the review. 
# ==============================================================================

# ---- 0. Load Libraries and Checkpoints --------------------------------------
library(Azimuth)
library(Seurat)
library(dplyr)
library(Azimuth)
library(patchwork)
library(ggplot2)

# Due to FastMNN requiring SCT assay split by sample to detect batches
# merged_all was splitted and is a Assay5 which PrepSCTFindMarkers() fails
# clean_SCT saved prior to SCTransform(), still SCT assay
merged_all <- readRDS("reproduce_figure5/data/merged_all_integrated.rds")
clean_SCT <- readRDS("reproduce_figure5/data/merged_all_sct.rds")

# ---- 1. Preprocessing --------------------------------------------------------
# copies over mnn and umap reduction to the SCTAsaay
clean_SCT[["mnn"]]  <- merged_all[["mnn"]]
clean_SCT[["umap"]] <- merged_all[["umap"]]
# copies over cluster assignment to the SCTAsaay
clean_SCT$seurat_clusters <- merged_all$seurat_clusters
Idents(clean_SCT) <- "seurat_clusters"

class(clean_SCT[["SCT"]])   # confirm it's still "SCTAssay"

# ---- 2. Find Markers ---------------------------------------------------------
# Corrects sequencing depth
clean_SCT <- PrepSCTFindMarkers(clean_SCT)

# Reports significant genes within each cluster
markers <- FindAllMarkers(clean_SCT, only.pos = TRUE)

# Find the top markers, group by cluster with 2 fold change (confident differential genes), only (+)
top_markers <- markers %>% 
  group_by(cluster) %>% 
  dplyr::filter(avg_log2FC > 1) %>% 
  dplyr::filter(p_val_adj < 0.05)

# Finds the top 10 marker gene per cluster
top10_markers <- top_markers %>% 
  group_by(cluster) %>% 
  slice_max(order_by = avg_log2FC, n = 10) %>%
  summarise(genes = paste(gene, collapse = ", "))

# ---- 3. Manual Annotation via CellMarker 3.0 ---------------------------------
# Using CellMarker 3.0
manual_labels <- c(
  "0" = "Myofibroblast",
  "1" = "Pericyte",
  "2" = "Smooth Muscle Cell",
  "3" = "Myeloid Dendritic Cell",
  "4" = "Alveolar Fibroblast",
  "5" = "Endothelial Cell",
  "6" = "Alveolar Epithelial Cell",
  "7" = "Natural Killer Cell",
  "8" = "Monocyte",
  "9" = "Myofibroblast",
  "10" = "Alveolar Macrophage",
  "11" = "Plasma Cell",
  "12" = "B Cell",
  "13" = "Mesothelial")

merged_all <- RenameIdents(merged_all, manual_labels)
merged_all$manual_labels <- Idents(merged_all)

manual_annot_plot <- DimPlot(merged_all, 
                             reduction = "umap", 
                             group.by = "manual_labels", 
                             label = TRUE, 
                             repel = TRUE) + NoLegend()
ggsave(
  filename = "reproduce_figure5/figures/manual_annotation_plot.png",
  plot = manual_annot_plot,
  width = 8,
  height = 6
)

# ---- 4. Automatic Annotation via Azimuth -------------------------------------
merged_all <- RunAzimuth(merged_all, reference = "lungref")

# SingleR annotation levels goes from broad to specific annotations
# tapply(X, INDEX, FUN) works through grouping, splits into X defined by INDEX and applies FUN
# Groups all Azimuth prediction by the cluster it belongs to 
# Function takes in the amount of label appears in the cluster (ex cluster 4 - Fibroblast 1000, Smooth muscle: 100), orders it and pulls the first label with the most amount
auto_cluster_majority <- tapply(merged_all$predicted.ann_level_3, merged_all$seurat_clusters,
                                function(x) names(sort(table(x), decreasing = TRUE))[1])
# Function takes the most frequent label and divide it by the number of cells and gives the % of which the cluster's cell got the same label.
auto_cluster_purity <- tapply(merged_all$predicted.ann_level_3, merged_all$seurat_clusters,
                              function(x) max(table(x)) / length(x))

auto_annot_plot <- DimPlot(merged_all, 
                           reduction = "umap", 
                           group.by = "predicted.ann_level_3", 
                           label = TRUE, 
                           repel = TRUE) + NoLegend()
ggsave(
  filename = "reproduce_figure5/figures/auto_annotation_plot.png",
  plot = auto_annot_plot,
  width = 8,
  height = 6
)

# ---- 5. Compare Manual and Automatic Annotations -----------------------------

lookup <- data.frame(
  cluster = names(auto_cluster_majority),
  azimuth_label = auto_cluster_majority,
  azimuth_purity = round(auto_cluster_purity, 2),
  manual_label = manual_labels[names(auto_cluster_majority)]
)
lookup

# ---- 5a. Test Cluster 9 ------------------------------------------------------
Idents(merged_all) <- "seurat_clusters"
markers_to_check <- c(
  "ACTA2", "TAGLN", "COL1A1", "POSTN", "CTHRC1", # myofibroblast / fibrotic markers
  "PLP1", "MPZ", "S100B" # neural / Schwann cell markers
)
cluster9_test_plot <- DotPlot(merged_all, features = markers_to_check) + RotatedAxis()
ggsave(
  filename = "reproduce_figure5/figures/cluster9_dotplot.png",
  plot = cluster9_test_plot,
  height = 8,
  width = 6
)
cluster9_plp1_plot <- FeaturePlot(merged_all, features = "PLP1", reduction = "umap")
ggsave(
  filename = "reproduce_figure5/figures/cluster9_plp1.png",
  plot = cluster9_plp1_plot,
  height = 8,
  width = 6
)
# ---- 6. Final Annotations ----------------------------------------------------
final_labels <- c(
  "0" = "Fibroblast",
  "1" = "Pericyte",
  "2" = "Smooth Muscle Cell",
  "3" = "Macrophage",
  "4" = "Alveolar Fibroblast",
  "5" = "Endothelial Cell",
  "6" = "Alveolar Epithelial Cell",
  "7" = "Natural Killer Cell",
  "8" = "Monocyte",
  "9" = "Myofibroblast",
  "10" = "Macrophage",
  "11" = "Plasma Cell",
  "12" = "B Cell",
  "13" = "Mesothelial")

# ---- 7. Rename Cluster using Final Labels ------------------------------------
Idents(merged_all) <- "seurat_clusters"
merged_all <- RenameIdents(merged_all, final_labels)
merged_all$cell_type <- Idents(merged_all)

# ---- 8. Visualizations -------------------------------------------------------

## ---- 8a. By Cell Type and By Condition --------------------------------------
annot_plot_celltype <- DimPlot(merged_all, reduction = "umap", group.by = "cell_type", label = TRUE, repel = TRUE) + NoLegend()
annot_plot_condition <- DimPlot(merged_all, reduction = "umap", group.by = "condition", repel = TRUE)

final_annot_plot <- (annot_plot_celltype + annot_plot_condition) + plot_annotation(title = "Final Annotated Cluster")
ggsave(
  filename = "reproduce_figure5/figures/final_annotation_plot.png",
  plot = final_annot_plot,
  width = 8,
  height = 6
)

## ---- 8b. Focus on COL1A1 ----------------------------------------------------
col1a1_plot <- FeaturePlot(merged_all, features = "COL1A1", reduction = "umap")
ggsave(
  filename = "reproduce_figure5/figures/col1a1_plot.png",
  plot = col1a1_plot,
  width = 8,
  height = 6
)

Idents(merged_all) <- "cell_type"   # make sure labels come from cell type, not seurat_clusters

col1a1_plot <- FeaturePlot(merged_all, features = "COL1A1", reduction = "umap",
                           order = TRUE, cols = c("lightgrey", "red"),
                           label = TRUE, repel = TRUE) & NoAxes()
ggsave(
  filename = "reproduce_figure5/figures/fig5b.png",
  plot = col1a1_plot,
  width = 8,
  height = 6
)
# Number of cells per condition
table(merged_all$condition)

## ---- 8c. # of Cell Types within each Condition ------------------------------
cell_props <- prop.table(table(merged_all$condition, merged_all$cell_type), margin = 1)
prop_df <- as.data.frame(cell_props)
colnames(prop_df) <- c("condition", "cell_type", "proportion")

celltype_per_condition <- ggplot(prop_df, aes(x = condition, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal()
ggsave(
  filename = "reproduce_figure5/figures/celltype_per_condition.png",
  plot = celltype_per_condition,
  height = 6,
  width = 8
)

# ---- 9. Checkpoint -----------------------------------------------------------
saveRDS(merged_all, file = "reproduce_figure5/data/merged_all_final_annotated.rds")
