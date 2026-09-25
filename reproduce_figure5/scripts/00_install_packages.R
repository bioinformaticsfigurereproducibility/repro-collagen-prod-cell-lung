# ============================================================================
# Install every package used across scripts/01-04 (+ 01b) in this project.
# Run once after cloning the repo, before running any numbered script.
# ============================================================================

cran_pkgs <- c(
  "Seurat",
  "sctransform",
  "ggplot2",
  "patchwork",
  "dplyr",
  "remotes",
  "here"
)

bioc_pkgs <- c(
  "SingleCellExperiment",
  "scran",
  "scater",
  "scuttle",
  "glmGamPoi",
  "batchelor"
)

github_pkgs <- c(
  "satijalab/seurat-wrappers",
  "satijalab/azimuth"
)

# --- CRAN ---
installed <- rownames(installed.packages())
new_cran <- setdiff(cran_pkgs, installed)
if (length(new_cran)) install.packages(new_cran)

# --- Bioconductor ---
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
new_bioc <- bioc_pkgs[!vapply(bioc_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(new_bioc)) BiocManager::install(new_bioc, update = FALSE, ask = FALSE)

# --- GitHub ---
for (pkg in github_pkgs) {
  pkg_name <- basename(pkg)
  if (!requireNamespace(pkg_name, quietly = TRUE)) {
    remotes::install_github(pkg, upgrade = "never")
  }
}

message("All packages installed.")

