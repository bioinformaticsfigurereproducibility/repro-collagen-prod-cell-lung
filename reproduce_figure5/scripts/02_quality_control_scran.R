# ==============================================================================
# scran-based QC, matching the paper's stated method:
#   "We excluded cells which were five median absolute deviations (MAD)
#    distant from the median value of library size, number of detected
#    genes, or mitochondrial gene proportion."
# ==============================================================================

# ---- 0. Load Packages --------------------------------------------------------
library(scran)
library(scater)
library(scuttle)
library(SingleCellExperiment)
library(Seurat)
library(here)

# ---- 1. Load raw (unfiltered) counts per sample ------------------------------
sample_dirs <- c(
  NML1_All = here("reproduce_figure5", "data", "NML1_All"), 
  NML2_All = here("reproduce_figure5", "data", "NML2_All"), 
  NML3_All = here("reproduce_figure5", "data", "NML3_All"),
  NML1_Lin = here("reproduce_figure5", "data", "NML1_Lin"), 
  NML2_Lin = here("reproduce_figure5", "data", "NML2_Lin"), 
  NML3_Lin = here("reproduce_figure5", "data", "NML3_Lin"),
  IPF1_All = here("reproduce_figure5", "data", "IPF1_All"), 
  IPF2_All = here("reproduce_figure5", "data", "IPF2_All"), 
  IPF3_All = here("reproduce_figure5", "data", "IPF3_All"),
  IPF1_Lin = here("reproduce_figure5", "data", "IPF1_Lin"), 
  IPF2_Lin = here("reproduce_figure5", "data", "IPF2_Lin"), 
  IPF3_Lin = here("reproduce_figure5", "data", "IPF3_Lin"),
  SCD1_All = here("reproduce_figure5", "data", "SCD1_All"), 
  SCD2_All = here("reproduce_figure5", "data", "SCD2_All"),
  SCD1_Lin = here("reproduce_figure5", "data", "SCD1_Lin"),
  SCD2_Lin = here("reproduce_figure5", "data", "SCD2_Lin")
)

counts_list <- lapply(sample_dirs, Read10X)

# ---- 2. Build one SingleCellExperiment, tagging each cell's sample -----------
# Apply a function of a list of more than 1
sce_list <- Map(function(counts, id) {
  colnames(counts) <- paste0(id, "_", colnames(counts)) # changes the column to the sample id and barcode
  SingleCellExperiment(assays = list(counts = counts)) # creates the sce object
}, counts_list, names(sample_dirs)) # input into function

sce <- do.call(cbind, sce_list) # combines the 16 sce objects into 1 sce object
sce$sample <- sub("_[ACGT]+.*$", "", colnames(sce)) # recover sample id per cell
sce$condition <- sub("[0-9].*$", "", sce$sample) # NML / IPF / SCD

# ---- 3. Per-cell QC metrics (library size, detected genes, mito%) ------------
is_mito <- grepl("^MT-", rownames(sce))
# computes the per-cell QC numbers
qc <- perCellQCMetrics(sce, subsets = list(Mito = is_mito))
colData(sce) <- cbind(colData(sce), qc) 

# ---- 4. 5-MAD outlier detection, per the paper's method --------------------
# paper excludes cells that were 5 median absolute deviation (MAD) distant from median
# type controls which directions counts as bad
# library size (right skewed)
qc_lib <- isOutlier(sce$sum, log = TRUE, type = "lower", nmads = 5, batch = sce$sample)
# number of genes (right skewed)
qc_genes <- isOutlier(sce$detected, log = TRUE,  type = "lower", nmads = 5, batch = sce$sample)
# mitochondrial gene (bounded to 0-100)
qc_mito <- isOutlier(sce$subsets_Mito_percent, log = FALSE, type = "higher", nmads = 5, batch = sce$sample)

discard <- qc_lib | qc_genes | qc_mito
sce$discard <- discard

# ---- 5. Summary of QC --------------------------------------------------------
message("QC SUMMARY (scran/scater, 5-MAD)")
message("Cells before filtering: ", ncol(sce)) # 83,704 cells
message("Cells flagged for removal: ", sum(discard),
        " (lib size: ", sum(qc_lib),
        ", detected genes: ", sum(qc_genes),
        ", mito%: ", sum(qc_mito), ")")

# ---- 6. Remove outliers ------------------------------------------------------
sce_filtered <- sce[, !discard] # removes cells that failed QC
message("Cells after filtering: ", ncol(sce_filtered)) # 81,474 cells

# ---- 7. Checkpoint -----------------------------------------------------------
saveRDS(sce_filtered, file = here("reproduce_figure5","data", "filtered_sce_scran_qc.rds"))
