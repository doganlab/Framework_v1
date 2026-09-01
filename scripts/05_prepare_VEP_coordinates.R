#!/usr/bin/env Rscript

# ============================================================
# Framework_v1
# 05_prepare_VEP_coordinates.R
#
# IMPORTANT BUILD NOTE
# --------------------
# Current Tambets et al. Nightingale mGWAS coordinates used by
# InstrumentEvidence are GRCh38 coordinates. VEP must therefore
# be run with a GRCh38 cache/assembly.
#
# This script prepares a deduplicated Ensembl-format VEP input
# from the candidate-linked InstrumentEvidence.tsv.
# ============================================================

options(stringsAsFactors = FALSE, warn = 1, scipen = 999)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
source(file.path(dirname(script_path), "..", "R", "bootstrap.R"))

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("The data.table package is required.", call. = FALSE)
}
library(data.table)

input_file <- INSTRUMENT_EVIDENCE_FILE
output_dir <- file.path(DIR_INTERMEDIATE, "variant_annotation")

output_file <- file.path(
  output_dir,
  "vep_coordinates_GRCh38.txt"
)

manifest_file <- file.path(
  output_dir,
  "VEPCoordinateManifest_GRCh38.tsv"
)

invalid_file <- file.path(
  output_dir,
  "InvalidCoordinateRecords_GRCh38.tsv"
)

conflict_file <- file.path(
  output_dir,
  "ConflictingVariantRecords_GRCh38.tsv"
)

qc_file <- file.path(
  output_dir,
  "VEP_coordinate_QC_GRCh38.tsv"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file, call. = FALSE)
}

dat <- fread(input_file, showProgress = FALSE)

required_columns <- c(
  "Metabolite",
  "StudyAccession",
  "SNP",
  "CHR",
  "BP",
  "EffectAllele",
  "OtherAllele"
)

missing_columns <- setdiff(required_columns, names(dat))

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

manifest <- dat[, .(
  Metabolite,
  StudyAccession,
  SNP = trimws(as.character(SNP)),
  CHR = sub(
    "^chr",
    "",
    trimws(as.character(CHR)),
    ignore.case = TRUE
  ),
  BP = suppressWarnings(as.integer(BP)),
  EffectAllele = toupper(trimws(as.character(EffectAllele))),
  OtherAllele = toupper(trimws(as.character(OtherAllele)))
)]

valid_chr <- manifest$CHR %in% c(
  as.character(1:22), "X", "Y", "MT", "M"
)

manifest[CHR == "M", CHR := "MT"]

valid_bp <- !is.na(manifest$BP) & manifest$BP > 0L
valid_snp <- !is.na(manifest$SNP) & nzchar(manifest$SNP)

valid_effect_allele <- (
  !is.na(manifest$EffectAllele) &
    grepl("^[ACGT]+$", manifest$EffectAllele)
)

valid_other_allele <- (
  !is.na(manifest$OtherAllele) &
    grepl("^[ACGT]+$", manifest$OtherAllele)
)

different_alleles <- (
  !is.na(manifest$EffectAllele) &
    !is.na(manifest$OtherAllele) &
    manifest$EffectAllele != manifest$OtherAllele
)

valid <- (
  valid_chr &
    valid_bp &
    valid_snp &
    valid_effect_allele &
    valid_other_allele &
    different_alleles
)

invalid <- manifest[!valid]
manifest_valid <- manifest[valid]

# One coordinate/allele record can support multiple metabolite exposures.
variant_manifest <- unique(
  manifest_valid[, .(
    SNP,
    CHR,
    BP,
    EffectAllele,
    OtherAllele
  )]
)

variant_manifest[, VariantKey := paste(
  CHR,
  BP,
  EffectAllele,
  OtherAllele,
  sep = ":"
)]

variant_key_counts <- variant_manifest[
  ,
  .(NVariantKeys = uniqueN(VariantKey)),
  by = SNP
]

conflicting_snps <- variant_key_counts[
  NVariantKeys > 1,
  SNP
]

conflicting <- variant_manifest[
  SNP %in% conflicting_snps
]

# Do not silently send ambiguous rsID-coordinate mappings to VEP.
# They are retained in a QC file and removed from the VEP manifest.
if (length(conflicting_snps) > 0) {
  warning(
    length(conflicting_snps),
    " rsID(s) have conflicting coordinate/allele records and will be excluded from VEP input."
  )
}

variant_manifest_clean <- variant_manifest[
  !SNP %in% conflicting_snps
]

variant_manifest_clean[, VariantKey := NULL]

# ------------------------------------------------------------
# Build sanity sentinel
# ------------------------------------------------------------
# rs964184 is a useful sentinel because the current mGWAS data
# carry it at chr11:116778201, consistent with GRCh38. If present
# at the GRCh37 position, stop rather than annotating against the
# wrong VEP cache.

if ("rs964184" %in% variant_manifest_clean$SNP) {
  sentinel <- variant_manifest_clean[SNP == "rs964184"][1]

  if (
    sentinel$CHR == "11" &&
      sentinel$BP == 116648917L
  ) {
    stop(
      paste0(
        "Genome-build sanity check failed: rs964184 is at ",
        "11:116648917, which corresponds to the older GRCh37 coordinate. ",
        "This VEP preparation script is configured for GRCh38."
      ),
      call. = FALSE
    )
  }

  if (
    sentinel$CHR == "11" &&
      sentinel$BP != 116778201L
  ) {
    warning(
      "rs964184 was found at an unexpected coordinate: ",
      sentinel$CHR, ":", sentinel$BP,
      ". Review genome build before running VEP."
    )
  }
}

# ------------------------------------------------------------
# Ensembl six-column input
# chromosome start end allele strand identifier
#
# For a single-base substitution, start=end=BP.
# Alleles are supplied as OTHER/EFFECT to preserve the exposure
# allele context while VEP determines the relevant allele.
# ------------------------------------------------------------

vep_input <- variant_manifest_clean[, .(
  CHR,
  START = BP,
  END = BP,
  ALLELE = paste0(OtherAllele, "/", EffectAllele),
  STRAND = "+",
  ID = SNP
)]

fwrite(
  vep_input,
  output_file,
  sep = "\t",
  col.names = FALSE,
  quote = FALSE,
  na = ""
)

fwrite(
  manifest_valid,
  manifest_file,
  sep = "\t",
  quote = FALSE,
  na = ""
)

fwrite(
  invalid,
  invalid_file,
  sep = "\t",
  quote = FALSE,
  na = ""
)

fwrite(
  conflicting,
  conflict_file,
  sep = "\t",
  quote = FALSE,
  na = ""
)

qc <- data.table(
  Metric = c(
    "InstrumentEvidence rows",
    "Valid instrument rows",
    "Invalid instrument rows",
    "Unique valid coordinate-allele records",
    "Unique rsIDs before conflict filtering",
    "Conflicting rsIDs excluded",
    "Unique VEP records",
    "Genome build for VEP"
  ),
  Value = c(
    as.character(nrow(dat)),
    as.character(nrow(manifest_valid)),
    as.character(nrow(invalid)),
    as.character(nrow(variant_manifest)),
    as.character(uniqueN(variant_manifest$SNP)),
    as.character(length(conflicting_snps)),
    as.character(nrow(vep_input)),
    "GRCh38"
  )
)

fwrite(qc, qc_file, sep = "\t", quote = FALSE)

cat("\nVEP coordinate input created.\n")
print(qc)

cat(
  "\nVEP input:      ", output_file,
  "\nManifest:       ", manifest_file,
  "\nInvalid:        ", invalid_file,
  "\nConflicts:      ", conflict_file,
  "\nQC:             ", qc_file,
  "\n\n",
  sep = ""
)

cat(
  "IMPORTANT: run VEP with --assembly GRCh38 and a matching GRCh38 cache.\n"
)
