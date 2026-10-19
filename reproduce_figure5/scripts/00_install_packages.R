# ============================================================================
# Install every package used across scripts/01-04 (+ 01b) in this project.
# Run once after cloning the repo, before running any numbered script.

# Run everything:Rscript 00_install_packages.R
# Run one section:Rscript 00_install_packages.R cran   (or bioc, github, check)
# ============================================================================

# ---- Helpers -----------------------------------------------------------------
# TRUE if the package is installed on disk.
# Checks the package folder instead of loading the package, so packages
# upgraded mid-install don't clash with versions already loaded in this session.
is_installed <- function(pkg) nzchar(system.file(package = pkg))

# Stop with a clear message if a package didn't install
require_installed <- function(pkg) {
  if (!is_installed(pkg)) {
    stop("Failed to install ", pkg, ". See the error messages above.")
  }
}

# TRUE only if the package is installed AND at the exact version
# (packageVersion() reads the DESCRIPTION file, does not load the package)
has_version <- function(pkg, version) {
  is_installed(pkg) && as.character(packageVersion(pkg)) == version
}

# ---- Installer tools ---------------------------------------------------------
if (!is_installed("remotes")) install.packages("remotes")
if (!is_installed("BiocManager")) install.packages("BiocManager")

# ---- Package lists (kept outside sections; the final check uses them) -------
cran_pkgs <- c(
  SeuratObject = "5.4.0",
  sctransform = "0.4.3",
  Seurat = "5.5.1",
  ggplot2 = "4.0.3",
  patchwork = "1.3.2",
  dplyr = "1.2.1",
  here = "1.0.2",
  R.utils = "2.13.0"
)

bioc_version <- "3.22"
bioc_pkgs <- c(
  SingleCellExperiment = "1.32.0",
  scuttle = "1.20.0",
  scran = "1.38.1",
  scater = "1.38.1",
  glmGamPoi = "1.22.0",
  batchelor = "1.26.0"
)

github_pkgs <- c(
  SeuratData  = "satijalab/seurat-data@3e51f44303069b64f5dc4d68e6a3d4a343f55c39",
  Banksy = "prabhakarlab/Banksy@5157d9cf8020ec49c99f18d833b143b69a496762",
  SeuratWrappers = "satijalab/seurat-wrappers@8df83435cf10ca76f28ddfc658ad3f6bec61824c",
  Azimuth = "satijalab/azimuth@ad5929686e005ae91d8452187a3a1d8c0563fdd5"
)

github_versions <- c(
  SeuratData = "0.2.2.9002",
  Banksy = "1.9.1",
  SeuratWrappers = "0.4.0",
  Azimuth = "0.5.1"
)

# ---- CRAN: exact versions ----------------------------------------------------
for (pkg in names(cran_pkgs)) {
  if (!has_version(pkg, cran_pkgs[[pkg]])) {
    remotes::install_version(pkg, version = cran_pkgs[[pkg]], upgrade = "never")
    require_installed(pkg)
  }
}

# ---- Bioconductor: pinned by release -----------------------------------------
new_bioc <- names(bioc_pkgs)[!vapply(names(bioc_pkgs), is_installed, logical(1))]
if (length(new_bioc)) {
  BiocManager::install(new_bioc, version = bioc_version,
                       update = FALSE, ask = FALSE)
}
for (pkg in names(bioc_pkgs)) require_installed(pkg)

# ---- GitHub: exact commits ---------------------------------------------------
for (pkg in names(github_pkgs)) {
  wanted_sha <- sub(".*@", "", github_pkgs[[pkg]])
  installed_sha <- if (is_installed(pkg)) {
    packageDescription(pkg)$RemoteSha
  } else NULL
  if (!identical(installed_sha, wanted_sha)) {
    remotes::install_github(github_pkgs[[pkg]], upgrade = "never")
    require_installed(pkg)
  }
}

# ---- Azimuth lung reference (used in 04_annotation.R) ------------------------
if (!is_installed("lungref.SeuratData")) {
  SeuratData::InstallData("lungref")
}
require_installed("lungref.SeuratData")

# ---- Check: does everything match? -------------------------------------------
expected <- c(cran_pkgs, bioc_pkgs, github_versions)
actual <- vapply(names(expected), function(p) {
  if (is_installed(p)) as.character(packageVersion(p)) else "NOT INSTALLED"
}, character(1))

if (any(actual == "NOT INSTALLED")) {
  print(data.frame(expected, actual)[actual == "NOT INSTALLED", ])
  stop("Some packages are missing (listed above).")
} else if (any(expected != actual)) {
  warning("Some versions differ from the ones used for the figures:")
  print(data.frame(expected, actual)[expected != actual, ])
} else {
  message("All packages installed at the expected versions.")
}