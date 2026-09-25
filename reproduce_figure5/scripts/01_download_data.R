# ==============================================================================
# Download the human lung scRNA-seq data (GEO: GSE132771) this project uses
# and arrange it into the data/<SAMPLE>_<Sort>/ layout Read10X() expects.
# ==============================================================================

library(here)

samples_name <- c(
  "NML1_Lin", "NML1_All", "NML2_Lin", "NML2_All", "NML3_Lin", "NML3_All",
  "IPF1_Lin", "IPF1_All", "IPF2_Lin", "IPF2_All", "IPF3_Lin", "IPF3_All",
  "SCD1_Lin", "SCD1_All", "SCD2_Lin", "SCD2_All"
)

raw_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE132nnn/GSE132771/suppl/GSE132771_RAW.tar"
raw_tar <- here("reproduce_figure5", "data", "GSE132771_RAW.tar")
extract_dir <- here("reproduce_figure5", "data", "_raw_extracted")

dir.create(here("reproduce_figure5", "data"), showWarnings = FALSE)
dir.create(extract_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(raw_tar)) {
  message("Downloading GSE132771_RAW.tar (~525 MB)")
  download.file(raw_url, destfile = raw_tar, mode = "wb", timeout = 1200)
}

untar(raw_tar, exdir = extract_dir)

all_files <- list.files(extract_dir, full.names = TRUE)

for (s in samples_name) {
  sample_dir <- here("reproduce_figure5", "data", s)
  dir.create(sample_dir, showWarnings = FALSE)
  
  matches <- all_files[grepl(paste0("_", s, "_"), all_files, fixed = FALSE)]
  
  if (length(matches) == 0) {
    warning("No files matched for sample '", s, "' - check reproduce_figure5/data/_raw_extracted manually.")
    next
  }
  
  for (f in matches) {
    target <- if (grepl("barcodes", f)) "barcodes.tsv.gz"
    else if (grepl("features|genes", f)) "features.tsv.gz"
    else if (grepl("matrix", f)) "matrix.mtx.gz"
    else NA
    if (!is.na(target)) file.copy(f, file.path(sample_dir, target), overwrite = TRUE)
  }
}

message("Completed downloading data.")