#!/usr/bin/env Rscript

# ============================================================
# Framework_v1
# 10_prepare_CPIC_benchmark.R
#
# Download CPIC pair_view from the official CPIC PostgREST API
# and prepare a compact, reproducible PGx benchmark table.
#
# IMPORTANT
# ---------
# This script DOES NOT automatically map CPIC individual drugs
# to Nature Medicine drug categories. That mapping should be
# reviewed explicitly before benchmark matching.
# ============================================================

options(stringsAsFactors = FALSE, warn = 1, scipen = 999)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
source(file.path(dirname(script_path), "..", "R", "bootstrap.R"))

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("The data.table package is required.", call. = FALSE)
}
library(data.table)

raw_dir <- file.path(DIR_RAW, "CPIC")
processed_dir <- DIR_PROCESSED
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

api_url <- "https://api.cpicpgx.org/v1/pair_view"

raw_file <- file.path(raw_dir, "CPIC_pair_view.csv")
benchmark_file <- file.path(processed_dir, "CPIC_PGx_Benchmark.tsv")
mapping_template_file <- file.path(
  processed_dir, "CPIC_to_NatMed_DrugCategory_Mapping.tsv"
)
qc_file <- file.path(processed_dir, "CPIC_PGx_Benchmark_QC.tsv")

cat("\n=== CPIC PGx benchmark preparation ===\n")
cat("Source: ", api_url, "\n", sep = "")

# Download as CSV. PostgREST supports text/csv via Accept header.
# download.file() cannot set headers portably on every R build,
# so curl is preferred when available.
curl_bin <- Sys.which("curl")

if (nzchar(curl_bin)) {
  status <- system2(
    curl_bin,
    args = c(
      "-L",
      "--fail",
      "--silent",
      "--show-error",
      "-H", shQuote("Accept: text/csv"),
      shQuote(api_url),
      "-o", shQuote(raw_file)
    )
  )
  if (!identical(status, 0L)) {
    stop("CPIC API download failed with curl.", call. = FALSE)
  }
} else {
  warning("curl executable not found; trying R download.file().")
  download.file(api_url, raw_file, mode = "wb", quiet = FALSE)
}

if (!file.exists(raw_file) || file.info(raw_file)$size == 0) {
  stop("Downloaded CPIC file is missing or empty.", call. = FALSE)
}

cpic <- fread(raw_file, na.strings = c("", "NA", "null"))

cat("Downloaded rows: ", nrow(cpic), "\n", sep = "")
cat("Columns: ", paste(names(cpic), collapse = ", "), "\n", sep = "")

# CPIC pair_view schema can evolve. Detect key columns defensively.
pick_col <- function(candidates, required = FALSE, label = NULL) {
  hit <- candidates[candidates %in% names(cpic)]
  if (length(hit) > 0) return(hit[1])
  if (required) {
    stop(
      "Could not identify required CPIC column",
      if (!is.null(label)) paste0(" for ", label) else "",
      ". Available columns: ", paste(names(cpic), collapse = ", "),
      call. = FALSE
    )
  }
  NA_character_
}

gene_col <- pick_col(
  c("genesymbol", "gene", "symbol"),
  required = TRUE,
  label = "gene symbol"
)
drug_col <- pick_col(
  c("drugname", "drug", "name"),
  required = TRUE,
  label = "drug name"
)
drugid_col <- pick_col(c("drugid", "drug_id"))
level_col <- pick_col(c("cpiclevel", "cpic_level"))
guideline_col <- pick_col(
  c("guidelinename", "guideline", "guideline_name")
)
guidelineid_col <- pick_col(c("guidelineid", "guideline_id"))
fda_col <- pick_col(
  c("pgxtesting", "fda_pgxtesting", "fdalabel", "fda_label")
)
clinpgx_col <- pick_col(
  c(
    "clinpgxlevel", "clinpgx_level",
    "pharmgkbca", "pharmgkb_level",
    "topclinicalannotationlevel"
  )
)

getv <- function(col) {
  if (is.na(col)) rep(NA_character_, nrow(cpic))
  else as.character(cpic[[col]])
}

benchmark <- data.table(
  Gene = toupper(trimws(getv(gene_col))),
  CPIC_Drug = trimws(getv(drug_col)),
  CPIC_DrugID = getv(drugid_col),
  CPIC_Level = getv(level_col),
  CPIC_Guideline = getv(guideline_col),
  CPIC_GuidelineID = getv(guidelineid_col),
  FDA_PGx_Field = getv(fda_col),
  ClinPGx_Field = getv(clinpgx_col)
)

benchmark[, CPIC_Drug_Normalized := tolower(trimws(CPIC_Drug))]
benchmark[, HasCPICGuideline :=
  (!is.na(CPIC_Guideline) & CPIC_Guideline != "") |
  (!is.na(CPIC_GuidelineID) & CPIC_GuidelineID != "")
]

# CPIC levels are retained exactly as supplied by CPIC.
# Do not equate all tracked pairs with guideline-backed pairs.
benchmark[, CPIC_A_or_B := CPIC_Level %in% c("A", "B")]

benchmark <- unique(benchmark)
setorder(benchmark, Gene, CPIC_Drug_Normalized)

fwrite(
  benchmark,
  benchmark_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# ------------------------------------------------------------
# Create a REVIEW TEMPLATE for NatMed category mapping.
# One row per CPIC drug. No automatic pharmacological mapping.
# ------------------------------------------------------------

mapping <- unique(
  benchmark[, .(
    CPIC_Drug,
    CPIC_Drug_Normalized,
    CPIC_DrugID
  )]
)

mapping[, NatMed_DrugCategory := NA_character_]
mapping[, NatMed_CategoryCode := NA_character_]
mapping[, MappingStatus := "UNREVIEWED"]
mapping[, MappingBasis := NA_character_]
mapping[, MappingNotes := NA_character_]

setorder(mapping, CPIC_Drug_Normalized)

fwrite(
  mapping,
  mapping_template_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

qc <- data.table(
  Metric = c(
    "CPIC pair_view rows downloaded",
    "Unique benchmark gene-drug pairs",
    "Unique CPIC genes",
    "Unique CPIC drugs",
    "CPIC level A pairs",
    "CPIC level B pairs",
    "CPIC level A or B pairs",
    "Pairs linked to a CPIC guideline",
    "Unique drugs requiring NatMed mapping review"
  ),
  Value = c(
    nrow(cpic),
    uniqueN(benchmark, by = c("Gene", "CPIC_Drug_Normalized")),
    uniqueN(benchmark$Gene),
    uniqueN(benchmark$CPIC_Drug_Normalized),
    benchmark[CPIC_Level == "A", .N],
    benchmark[CPIC_Level == "B", .N],
    benchmark[CPIC_A_or_B == TRUE, .N],
    benchmark[HasCPICGuideline == TRUE, .N],
    nrow(mapping)
  )
)

fwrite(qc, qc_file, sep = "\t", quote = FALSE)

cat("\n============================================================\n")
cat("CPIC benchmark preparation completed\n")
cat("============================================================\n\n")
print(qc)

cat(
  "\nRaw CPIC snapshot:\n  ", raw_file,
  "\n\nProcessed benchmark:\n  ", benchmark_file,
  "\n\nDrug-category mapping template:\n  ", mapping_template_file,
  "\n\nQC:\n  ", qc_file,
  "\n\n",
  sep = ""
)

cat(
  "NEXT STEP:\n",
  "Review CPIC_to_NatMed_DrugCategory_Mapping.tsv before any\n",
  "gene-drug recovery analysis. Do not automatically interpret\n",
  "all CPIC-tracked pairs as guideline-backed positive controls.\n"
)
