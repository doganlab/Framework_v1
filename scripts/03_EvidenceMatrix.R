#!/usr/bin/env Rscript

# ============================================================
# Framework_v1
# 03_EvidenceMatrix.R
#
# Nature Medicine version
#
# Purpose
# -------
# Freeze the pharmacometabolomic evidence layer and define the
# primary metabolite-drug candidate set WITHOUT a composite score.
#
# Inputs
# ------
# PharmacoMetabolomicEvidence.tsv
#
# Principal outputs
# -----------------
# EvidenceMatrix.tsv                    : all 8,874 evidence rows
# DrugClassCandidateNetwork.tsv         : primary evidence rows
# Tier1MetaboliteDrugEvidence.tsv       : robust primary rows
# Tier2MetaboliteDrugEvidence.tsv       : high-confidence primary rows
# DirectionallyDiscordantEvidence.tsv   : significant but opposite interpretation
# EvidenceMatrix_summary.tsv            : QC / manuscript-ready counts
# CandidateDrugCategorySummary.tsv
# CandidateMetaboliteSummary.tsv
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
# Paths
# ------------------------------------------------------------

input_file <- PHARMACOMETABOLOMICS_FILE

tier1_file <- file.path(
  DIR_PROCESSED,
  "Tier1MetaboliteDrugEvidence.tsv"
)

tier2_file <- file.path(
  DIR_PROCESSED,
  "Tier2MetaboliteDrugEvidence.tsv"
)

discordant_file <- file.path(
  DIR_PROCESSED,
  "DirectionallyDiscordantEvidence.tsv"
)

summary_file <- file.path(
  DIR_PROCESSED,
  "EvidenceMatrix_summary.tsv"
)

drug_summary_file <- file.path(
  DIR_PROCESSED,
  "CandidateDrugCategorySummary.tsv"
)

metabolite_summary_file <- file.path(
  DIR_PROCESSED,
  "CandidateMetaboliteSummary.tsv"
)

create_directory(DIR_PROCESSED)
check_file_exists(input_file)

# ------------------------------------------------------------
# Read and validate evidence
# ------------------------------------------------------------

cat("\n=== Building score-free evidence matrix ===\n")

evidence <- fread(input_file, showProgress = FALSE)

required_cols <- c(
  "Metabolite",
  "DrugCategory",
  "CategoryCode",
  "MR_beta",
  "MR_p",
  "MR_FDR",
  "PP_H4",
  "M2_Beta",
  "M2_P",
  "M2Significant",
  "RiskOpposingDirection",
  "M1M2Robust",
  "PrimaryPharmacoGeneticEvidence",
  "EvidenceTier",
  "CardiovascularATCGroup"
)

check_columns(evidence, required_cols)

# Each Nature Medicine trait should have one row for each of 87 drug categories.
n_metabolites <- uniqueN(evidence$Metabolite)
n_drugs <- uniqueN(evidence$DrugCategory)
expected_rows <- n_metabolites * n_drugs

if (nrow(evidence) != expected_rows) {
  stop_framework(
    paste0(
      "Evidence matrix is not a complete metabolite x drug-category grid. ",
      "Observed rows: ", nrow(evidence),
      "; expected: ", expected_rows,
      " (", n_metabolites, " x ", n_drugs, ")."
    )
  )
}

duplicate_edges <- evidence[
  ,
  .N,
  by = .(Metabolite, DrugCategory)
][N > 1]

if (nrow(duplicate_edges) > 0) {
  stop_framework(
    paste0(
      "Duplicate metabolite-drug edges detected: ",
      nrow(duplicate_edges)
    )
  )
}

# ------------------------------------------------------------
# Full evidence matrix
#
# No FinalPriorityScore is calculated. Independent evidence
# columns and pre-defined evidence tiers are preserved.
# ------------------------------------------------------------

evidence[, CandidateSelected :=
  PrimaryPharmacoGeneticEvidence %in% TRUE
]

evidence[, CandidateSelectionReason := fifelse(
  EvidenceTier == "Tier1_Robust",
  "PrimaryEvidence;Tier1_Robust",
  fifelse(
    EvidenceTier == "Tier2_HighConfidence",
    "PrimaryEvidence;Tier2_HighConfidence",
    fifelse(
      EvidenceTier == "DirectionallyDiscordant",
      "SignificantEvidence;DirectionallyDiscordant",
      "NotPrimary"
    )
  )
)]

setorder(
  evidence,
  -CandidateSelected,
  EvidenceTier,
  MR_FDR,
  M2_P,
  Metabolite,
  DrugCategory
)

write_tsv_safe(
  evidence,
  EVIDENCE_MATRIX_FILE
)

# ------------------------------------------------------------
# Primary candidate network
# ------------------------------------------------------------

candidate_network <- evidence[
  PrimaryPharmacoGeneticEvidence %in% TRUE
]

if (nrow(candidate_network) == 0) {
  stop_framework("Primary candidate network contains zero rows.")
}

if (any(candidate_network$EvidenceTier %in%
        c("DirectionallyDiscordant", "NotPrimary"))) {
  stop_framework(
    "Primary candidate network contains non-primary evidence tiers."
  )
}

setorder(
  candidate_network,
  EvidenceTier,
  MR_FDR,
  M2_P,
  Metabolite,
  DrugCategory
)

write_tsv_safe(
  candidate_network,
  DRUG_CANDIDATE_NETWORK_FILE
)

tier1 <- candidate_network[EvidenceTier == "Tier1_Robust"]
tier2 <- candidate_network[EvidenceTier == "Tier2_HighConfidence"]
discordant <- evidence[EvidenceTier == "DirectionallyDiscordant"]

write_tsv_safe(tier1, tier1_file)
write_tsv_safe(tier2, tier2_file)
write_tsv_safe(discordant, discordant_file)

# ------------------------------------------------------------
# Manuscript / QC summaries
# ------------------------------------------------------------

summary <- data.table(
  Metric = c(
    "All metabolite-drug edges",
    "Metabolites in full matrix",
    "Drug categories in full matrix",
    "M2-significant edges",
    "M2-significant and risk-opposing edges",
    "Primary pharmacogenetic evidence edges",
    "Primary metabolites",
    "Primary drug categories",
    "Tier1 robust edges",
    "Tier2 high-confidence edges",
    "Directionally discordant edges",
    "Primary cardiovascular-ATC edges",
    "Primary cardiovascular-ATC drug categories"
  ),
  Value = c(
    nrow(evidence),
    uniqueN(evidence$Metabolite),
    uniqueN(evidence$DrugCategory),
    evidence[M2Significant %in% TRUE, .N],
    evidence[
      M2Significant %in% TRUE &
        RiskOpposingDirection %in% TRUE,
      .N
    ],
    nrow(candidate_network),
    uniqueN(candidate_network$Metabolite),
    uniqueN(candidate_network$DrugCategory),
    nrow(tier1),
    nrow(tier2),
    nrow(discordant),
    candidate_network[
      CardiovascularATCGroup %in% TRUE,
      .N
    ],
    uniqueN(
      candidate_network[
        CardiovascularATCGroup %in% TRUE,
        DrugCategory
      ]
    )
  )
)

write_tsv_safe(summary, summary_file)

drug_summary <- candidate_network[
  ,
  .(
    PrimaryEdges = .N,
    Tier1Edges = sum(EvidenceTier == "Tier1_Robust"),
    Tier2Edges = sum(EvidenceTier == "Tier2_HighConfidence"),
    UniqueMetabolites = uniqueN(Metabolite),
    CardiovascularATCGroup = any(
      CardiovascularATCGroup %in% TRUE,
      na.rm = TRUE
    ),
    MedianAbsMRBeta = median(abs(MR_beta), na.rm = TRUE),
    MedianAbsM2Beta = median(abs(M2_Beta), na.rm = TRUE)
  ),
  by = .(
    DrugCategory,
    CategoryCode,
    DrugName,
    AnatomicalMainGroup,
    GroupLetter
  )
]

setorder(drug_summary, -PrimaryEdges, -Tier1Edges, DrugCategory)
write_tsv_safe(drug_summary, drug_summary_file)

metabolite_summary <- candidate_network[
  ,
  .(
    PrimaryEdges = .N,
    Tier1Edges = sum(EvidenceTier == "Tier1_Robust"),
    Tier2Edges = sum(EvidenceTier == "Tier2_HighConfidence"),
    UniqueDrugCategories = uniqueN(DrugCategory),
    CardiovascularDrugCategories = uniqueN(
      DrugCategory[CardiovascularATCGroup %in% TRUE]
    ),
    MR_beta = first(MR_beta),
    MR_p = first(MR_p),
    MR_FDR = first(MR_FDR),
    PP_H4 = first(PP_H4),
    MR_Direction = first(MR_Direction)
  ),
  by = .(
    Metabolite,
    NatMedAbbreviation,
    NatMedExplanation,
    MetaboliteGroup
  )
]

setorder(metabolite_summary, -PrimaryEdges, MR_FDR, Metabolite)
write_tsv_safe(metabolite_summary, metabolite_summary_file)

cat("\nEvidence matrix completed.\n")
print(summary)

cat(
  "\nFull matrix:      ", EVIDENCE_MATRIX_FILE,
  "\nPrimary network:  ", DRUG_CANDIDATE_NETWORK_FILE,
  "\nTier 1:           ", tier1_file,
  "\nTier 2:           ", tier2_file,
  "\nDiscordant:       ", discordant_file,
  "\nDrug summary:     ", drug_summary_file,
  "\nMetabolite summary:", metabolite_summary_file,
  "\n\n",
  sep = ""
)
