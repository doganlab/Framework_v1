#!/usr/bin/env Rscript

# ============================================================
# Framework_v1
# 04_VariantEvidence.R
#
# Nature Medicine version
#
# Purpose
# -------
# Propagate the primary metabolite-drug evidence to the
# independent harmonized MR instruments:
#
# SNP -> metabolite -> drug category
#
# This script does NOT infer a direct SNP-drug causal effect.
# Each row is an evidence chain linking a metabolite instrument
# to a metabolite-drug category association.
#
# Inputs
# ------
# DrugClassCandidateNetwork.tsv
# data/intermediate/MR/harmonized/*_harmonized_for_MR_rsID.tsv
#
# Outputs
# -------
# InstrumentEvidence.tsv
# MetaboliteInstrumentSummary.tsv
# VariantDrugClassEvidence.tsv
# EvidenceMatrix_withVariants.tsv
# DrugClassCandidateNetwork_withVariants.tsv
# VariantEvidence_summary.tsv
#
# InstrumentEvidence.tsv is intentionally restricted to
# metabolites represented in the primary candidate network.
# This keeps the subsequent VEP annotation focused on variants
# participating in the final evidence chains.
# ============================================================

options(stringsAsFactors = FALSE, warn = 1, scipen = 999)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)

source(file.path(dirname(script_path), "..", "R", "bootstrap.R"))
source(file.path(PROJECT_ROOT, "R", "io.R"))
source(file.path(PROJECT_ROOT, "R", "utils.R"))

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("The data.table package is required.", call. = FALSE)
}

library(data.table)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

read_one_harmonized_file <- function(file) {

  x <- fread(file, showProgress = FALSE)

  required <- c(
    "SNP",
    "beta.exposure",
    "se.exposure",
    "effect_allele.exposure",
    "other_allele.exposure",
    "eaf.exposure",
    "pval.exposure",
    "beta.outcome",
    "se.outcome",
    "pval.outcome",
    "CHR_exposure",
    "BP_exposure",
    "harmonisation",
    "palindromic"
  )

  check_columns(x, required)

  metabolite <- sub(
    HARMONIZED_INSTRUMENT_PATTERN,
    "",
    basename(file)
  )

  x[, Metabolite := metabolite]

  x
}

standardize_instrument_evidence <- function(
    harmonized,
    metabolite_metadata
) {

  # One exposure-specific record per metabolite-SNP.
  setorder(harmonized, Metabolite, SNP, pval.exposure)

  harmonized <- unique(
    harmonized,
    by = c("Metabolite", "SNP")
  )

  out <- harmonized[, .(
    Metabolite,
    SNP = as.character(SNP),
    CHR = as.character(CHR_exposure),
    BP = as.integer(BP_exposure),
    EffectAllele = toupper(as.character(effect_allele.exposure)),
    OtherAllele = toupper(as.character(other_allele.exposure)),
    EffectAlleleFrequency = as.numeric(eaf.exposure),
    BetaExposure = as.numeric(beta.exposure),
    SEExposure = as.numeric(se.exposure),
    PExposure = as.numeric(pval.exposure),
    BetaOutcome = as.numeric(beta.outcome),
    SEOutcome = as.numeric(se.outcome),
    POutcome = as.numeric(pval.outcome),
    HarmonizationStatus = as.character(harmonisation),
    Palindromic = as.logical(palindromic)
  )]

  out[, FStatistic := (BetaExposure / SEExposure)^2]

  out[, WeakInstrument :=
    is.na(FStatistic) | FStatistic <= 10
  ]

  out[, InstrumentStrengthClass := fifelse(
    is.na(FStatistic),
    "Unknown",
    fifelse(
      FStatistic > 10,
      "Strong",
      "Weak"
    )
  )]

  out <- merge(
    out,
    metabolite_metadata,
    by = "Metabolite",
    all.x = TRUE,
    sort = FALSE
  )

  if (any(is.na(out$StudyAccession) | out$StudyAccession == "")) {
    stop_framework(
      "At least one instrument failed to map to a metabolite study accession."
    )
  }

  preferred <- c(
    "Metabolite",
    "StudyAccession",
    "NatMedAbbreviation",
    "NatMedExplanation",
    "MetaboliteGroup",
    "SNP",
    "CHR",
    "BP",
    "EffectAllele",
    "OtherAllele",
    "EffectAlleleFrequency",
    "BetaExposure",
    "SEExposure",
    "PExposure",
    "BetaOutcome",
    "SEOutcome",
    "POutcome",
    "FStatistic",
    "InstrumentStrengthClass",
    "WeakInstrument",
    "HarmonizationStatus",
    "Palindromic"
  )

  keep <- intersect(preferred, names(out))
  out <- out[, ..keep]

  setorder(out, Metabolite, PExposure, SNP)

  out[]
}

# ------------------------------------------------------------
# Validate input paths
# ------------------------------------------------------------

check_file_exists(DRUG_CANDIDATE_NETWORK_FILE)
check_dir_exists(HARMONIZED_INSTRUMENTS_DIR)

create_directory(DIR_PROCESSED)

cat("\n=== Building variant-metabolite-drug evidence ===\n")

candidate_network <- fread(
  DRUG_CANDIDATE_NETWORK_FILE,
  showProgress = FALSE
)

required_candidate_cols <- c(
  "Metabolite",
  "study_accession",
  "NatMedAbbreviation",
  "MR_beta",
  "MR_p",
  "MR_FDR",
  "PP_H4",
  "DrugCategory",
  "CategoryCode",
  "M2_Beta",
  "M2_P",
  "M2Significant",
  "RiskOpposingDirection",
  "M1M2Robust",
  "EvidenceTier",
  "PrimaryPharmacoGeneticEvidence"
)

check_columns(candidate_network, required_candidate_cols)

if (any(!candidate_network$PrimaryPharmacoGeneticEvidence %in% TRUE)) {
  stop_framework(
    "DrugClassCandidateNetwork contains non-primary rows."
  )
}

candidate_metabolites <- sort(
  unique(candidate_network$Metabolite)
)

cat("Primary edges:       ", nrow(candidate_network), "\n", sep = "")
cat("Primary metabolites: ", length(candidate_metabolites), "\n", sep = "")
cat(
  "Primary drug classes:",
  uniqueN(candidate_network$DrugCategory),
  "\n"
)

# ------------------------------------------------------------
# Read only candidate-linked harmonized instrument files
# ------------------------------------------------------------

files <- list.files(
  HARMONIZED_INSTRUMENTS_DIR,
  pattern = HARMONIZED_INSTRUMENT_PATTERN,
  full.names = TRUE
)

if (length(files) == 0) {
  stop_framework(
    paste0(
      "No harmonized instrument files found in: ",
      HARMONIZED_INSTRUMENTS_DIR
    )
  )
}

file_metabolites <- sub(
  HARMONIZED_INSTRUMENT_PATTERN,
  "",
  basename(files)
)

keep_files <- file_metabolites %in% candidate_metabolites
files <- files[keep_files]
file_metabolites <- file_metabolites[keep_files]

missing_metabolites <- setdiff(
  candidate_metabolites,
  file_metabolites
)

if (length(missing_metabolites) > 0) {
  stop_framework(
    paste0(
      "Missing harmonized instrument files for candidate metabolites:\n",
      paste(missing_metabolites, collapse = "\n")
    )
  )
}

harmonized_list <- lapply(
  files,
  read_one_harmonized_file
)

harmonized <- rbindlist(
  harmonized_list,
  fill = TRUE,
  use.names = TRUE
)

cat(
  "Harmonized instrument rows read: ",
  format(nrow(harmonized), big.mark = ","),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# Metabolite metadata
# ------------------------------------------------------------

metabolite_metadata <- unique(
  candidate_network[, .(
    Metabolite,
    StudyAccession = study_accession,
    NatMedAbbreviation,
    NatMedExplanation,
    MetaboliteGroup
  )],
  by = "Metabolite"
)

if (nrow(metabolite_metadata) != length(candidate_metabolites)) {
  stop_framework(
    "Candidate metabolite metadata is not one-to-one."
  )
}

# ------------------------------------------------------------
# InstrumentEvidence
# ------------------------------------------------------------

instrument_evidence <- standardize_instrument_evidence(
  harmonized = harmonized,
  metabolite_metadata = metabolite_metadata
)

if (anyDuplicated(
  instrument_evidence[, paste(Metabolite, SNP, sep = "::")]
)) {
  stop_framework(
    "Duplicate metabolite-SNP pairs remain in InstrumentEvidence."
  )
}

write_tsv_safe(
  instrument_evidence,
  INSTRUMENT_EVIDENCE_FILE
)

instrument_summary <- instrument_evidence[
  ,
  .(
    NInstruments = .N,
    NUniqueSNPs = uniqueN(SNP),
    NStrongInstruments = sum(
      InstrumentStrengthClass == "Strong",
      na.rm = TRUE
    ),
    NWeakInstruments = sum(
      InstrumentStrengthClass == "Weak",
      na.rm = TRUE
    ),
    MedianF = median(FStatistic, na.rm = TRUE),
    MinF = min(FStatistic, na.rm = TRUE),
    MaxF = max(FStatistic, na.rm = TRUE),
    NPalindromic = sum(Palindromic %in% TRUE, na.rm = TRUE)
  ),
  by = .(
    Metabolite,
    StudyAccession,
    NatMedAbbreviation
  )
]

setorder(instrument_summary, -NInstruments, Metabolite)

write_tsv_safe(
  instrument_summary,
  METABOLITE_INSTRUMENT_SUMMARY_FILE
)

# ------------------------------------------------------------
# Attach instrument summaries to metabolite-drug evidence
# ------------------------------------------------------------

evidence_matrix <- fread(
  EVIDENCE_MATRIX_FILE,
  showProgress = FALSE
)

evidence_matrix_with_variants <- merge(
  evidence_matrix,
  instrument_summary[
    ,
    .(
      Metabolite,
      NInstruments,
      NStrongInstruments,
      NWeakInstruments,
      MedianF,
      MinF,
      MaxF
    )
  ],
  by = "Metabolite",
  all.x = TRUE,
  sort = FALSE
)

candidate_network_with_variants <- merge(
  candidate_network,
  instrument_summary[
    ,
    .(
      Metabolite,
      NInstruments,
      NStrongInstruments,
      NWeakInstruments,
      MedianF,
      MinF,
      MaxF
    )
  ],
  by = "Metabolite",
  all.x = TRUE,
  sort = FALSE
)

write_tsv_safe(
  evidence_matrix_with_variants,
  EVIDENCE_MATRIX_WITH_VARIANTS_FILE
)

write_tsv_safe(
  candidate_network_with_variants,
  CANDIDATE_NETWORK_WITH_VARIANTS_FILE
)

# ------------------------------------------------------------
# SNP -> metabolite -> drug-category propagation
# ------------------------------------------------------------

variant_drug <- merge(
  instrument_evidence,
  candidate_network,
  by = "Metabolite",
  allow.cartesian = TRUE,
  all = FALSE,
  sort = FALSE
)

# Avoid duplicated metadata fields from both sides.
# Candidate-network StudyAccession is study_accession, while
# InstrumentEvidence carries the canonical StudyAccession.
drop_if_present <- intersect(
  c(
    "study_accession",
    "NatMedAbbreviation.y",
    "NatMedExplanation.y",
    "MetaboliteGroup.y"
  ),
  names(variant_drug)
)

if (length(drop_if_present) > 0) {
  variant_drug[, (drop_if_present) := NULL]
}

rename_if_present <- function(dt, old, new) {
  if (old %in% names(dt) && !(new %in% names(dt))) {
    setnames(dt, old, new)
  }
}

rename_if_present(
  variant_drug,
  "NatMedAbbreviation.x",
  "NatMedAbbreviation"
)

rename_if_present(
  variant_drug,
  "NatMedExplanation.x",
  "NatMedExplanation"
)

rename_if_present(
  variant_drug,
  "MetaboliteGroup.x",
  "MetaboliteGroup"
)

# ------------------------------------------------------------
# Expansion validation
# ------------------------------------------------------------

instrument_counts <- instrument_evidence[
  ,
  .N,
  by = Metabolite
]

edge_counts <- candidate_network[
  ,
  .N,
  by = Metabolite
]

expected_expansion <- merge(
  instrument_counts,
  edge_counts,
  by = "Metabolite",
  suffixes = c("_instrument", "_edge")
)

expected_rows <- sum(
  expected_expansion$N_instrument *
    expected_expansion$N_edge
)

if (nrow(variant_drug) != expected_rows) {
  stop_framework(
    paste0(
      "Variant evidence expansion mismatch. Observed rows: ",
      nrow(variant_drug),
      "; expected: ",
      expected_rows,
      "."
    )
  )
}

required_variant_cols <- c(
  "SNP",
  "CHR",
  "BP",
  "EffectAllele",
  "OtherAllele",
  "Metabolite",
  "DrugCategory",
  "MR_beta",
  "MR_FDR",
  "PP_H4",
  "M2_Beta",
  "M2_P",
  "EvidenceTier"
)

check_columns(variant_drug, required_variant_cols)

setorder(
  variant_drug,
  EvidenceTier,
  MR_FDR,
  M2_P,
  Metabolite,
  DrugCategory,
  PExposure,
  SNP
)

write_tsv_safe(
  variant_drug,
  VARIANT_DRUG_EVIDENCE_FILE
)

# ------------------------------------------------------------
# QC summary
# ------------------------------------------------------------

summary_file <- file.path(
  DIR_PROCESSED,
  "VariantEvidence_summary.tsv"
)

summary <- data.table(
  Metric = c(
    "Primary metabolite-drug edges",
    "Candidate-linked metabolites",
    "Candidate-linked drug categories",
    "Unique metabolite-SNP instrument pairs",
    "Unique SNP rsIDs",
    "Strong instrument rows",
    "Weak instrument rows",
    "Variant-metabolite-drug evidence rows",
    "Tier1 variant evidence rows",
    "Tier2 variant evidence rows",
    "Cardiovascular-ATC variant evidence rows"
  ),
  Value = c(
    nrow(candidate_network),
    uniqueN(candidate_network$Metabolite),
    uniqueN(candidate_network$DrugCategory),
    nrow(instrument_evidence),
    uniqueN(instrument_evidence$SNP),
    instrument_evidence[
      InstrumentStrengthClass == "Strong",
      .N
    ],
    instrument_evidence[
      InstrumentStrengthClass == "Weak",
      .N
    ],
    nrow(variant_drug),
    variant_drug[EvidenceTier == "Tier1_Robust", .N],
    variant_drug[EvidenceTier == "Tier2_HighConfidence", .N],
    variant_drug[CardiovascularATCGroup %in% TRUE, .N]
  )
)

write_tsv_safe(summary, summary_file)

cat("\nVariant evidence completed.\n")
print(summary)

cat(
  "\nInstrument evidence: ",
  INSTRUMENT_EVIDENCE_FILE,
  "\nVariant-drug evidence: ",
  VARIANT_DRUG_EVIDENCE_FILE,
  "\nSummary: ",
  summary_file,
  "\n\n",
  sep = ""
)

cat(
  "Next step: run 05_prepare_VEP_coordinates.R using ",
  "InstrumentEvidence.tsv.\n",
  sep = ""
)
