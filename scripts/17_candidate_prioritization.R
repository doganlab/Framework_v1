#!/usr/bin/env Rscript

# ============================================================
# Framework_v1
# 21_candidate_prioritization.R
#
# Purpose:
#   Separate established benchmark recovery from framework-
#   prioritized candidate gene-drug modules.
#
# Design principle:
#   This script does NOT create an arbitrary composite score.
#   It applies transparent evidence gates and preserves every
#   evidence dimension for downstream interpretation.
#
# Candidate classes
# -----------------
# Established benchmark:
#   CPIC Level A/B OR ChEMBL single-protein direct MoA.
#
# Candidate Tier A:
#   cardiovascular drug category
#   + >=2 Tier-1 SNPs
#   + >=5 Tier-1 metabolites
#   + functional variant support
#   + DGIdb support
#   + NOT established benchmark
#
# Candidate Tier B:
#   same internal convergence and functional support,
#   but without requiring DGIdb support
#   + NOT established benchmark.
#
# Candidate Tier C:
#   >=2 Tier-1 SNPs + >=5 Tier-1 metabolites
#   + NOT established benchmark.
#
# IMPORTANT:
#   "Candidate" means hypothesis-generating framework output.
#   It does NOT mean novel, causal, clinically actionable,
#   or therapeutically beneficial.
# ============================================================

options(stringsAsFactors = FALSE, warn = 1, scipen = 999)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
source(file.path(dirname(script_path), "..", "R", "bootstrap.R"))

for (pkg in c("data.table")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required R package missing: ", pkg, call. = FALSE)
  }
}
library(data.table)

INPUT <- file.path(
  PROJECT_ROOT, "results", "integrated_evidence",
  "IntegratedGeneDrugEvidence.tsv"
)
OUT_DIR <- file.path(PROJECT_ROOT, "results", "candidate_prioritization")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(INPUT)) stop("Missing: ", INPUT, call. = FALSE)

x <- fread(INPUT, sep = "\t", na.strings = c("", "NA"), showProgress = FALSE)

required <- c(
  "Gene", "NatMed_DrugCategory", "NatMed_CategoryCode",
  "CardiovascularATCGroup",
  "UniqueSNPs", "UniqueMetabolites",
  "Tier1SNPs", "Tier1Metabolites",
  "Tier2SNPs", "Tier2Metabolites",
  "CodingSNPs", "SpliceRelatedSNPs",
  "HighImpactSNPs", "ModerateImpactSNPs",
  "OverlappingGeneSNPs", "NearestGeneSNPs",
  "MinMRFDR", "MaxPPH4", "MinM2P",
  "SNPs", "Metabolites",
  "FunctionalVariantInternalSupport",
  "CPICTrackedPair", "CPICLevelABSupport",
  "DGIdbSupport",
  "ChEMBLMoASupport", "ChEMBLSingleProteinMoASupport",
  "ExternalDatabasePattern"
)
miss <- setdiff(required, names(x))
if (length(miss)) {
  stop("IntegratedGeneDrugEvidence missing: ",
       paste(miss, collapse = ", "), call. = FALSE)
}

as_flag <- function(z) {
  if (is.logical(z)) return(!is.na(z) & z)
  toupper(trimws(as.character(z))) %in% c("TRUE", "YES", "Y", "1")
}

flag_cols <- c(
  "CardiovascularATCGroup",
  "FunctionalVariantInternalSupport",
  "CPICTrackedPair", "CPICLevelABSupport",
  "DGIdbSupport", "ChEMBLMoASupport",
  "ChEMBLSingleProteinMoASupport"
)
for (cc in flag_cols) set(x, j = cc, value = as_flag(x[[cc]]))

# ------------------------------------------------------------
# Evidence gates
# ------------------------------------------------------------

x[, EstablishedBenchmark :=
    CPICLevelABSupport | ChEMBLSingleProteinMoASupport]

x[, InternalConvergence :=
    CardiovascularATCGroup &
    Tier1SNPs >= 2L &
    Tier1Metabolites >= 5L]

x[, CandidateTierA :=
    InternalConvergence &
    FunctionalVariantInternalSupport &
    DGIdbSupport &
    !EstablishedBenchmark]

x[, CandidateTierB :=
    InternalConvergence &
    FunctionalVariantInternalSupport &
    !EstablishedBenchmark]

x[, CandidateTierC :=
    InternalConvergence &
    !EstablishedBenchmark]

x[, CandidateClass := fifelse(
  EstablishedBenchmark, "Established benchmark",
  fifelse(
    CandidateTierA, "Tier A: convergent + functional + DGIdb",
    fifelse(
      CandidateTierB, "Tier B: convergent + functional",
      fifelse(
        CandidateTierC, "Tier C: convergent internal evidence",
        "Not prioritized"
      )
    )
  )
)]

# Descriptive breadth only — explicitly NOT a ranking score.
x[, InternalEvidenceBreadth := rowSums(cbind(
  Tier1SNPs >= 2L,
  Tier1Metabolites >= 5L,
  FunctionalVariantInternalSupport,
  CodingSNPs > 0L,
  SpliceRelatedSNPs > 0L,
  HighImpactSNPs > 0L,
  ModerateImpactSNPs > 0L
), na.rm = TRUE)]

# ------------------------------------------------------------
# Transparent ordering (not scoring)
# ------------------------------------------------------------

order_candidates <- function(dt) {
  z <- copy(dt)
  setorder(
    z,
    -Tier1Metabolites,
    -Tier1SNPs,
    -UniqueMetabolites,
    -UniqueSNPs,
    Gene,
    NatMed_DrugCategory
  )
  z
}

tierA <- order_candidates(x[CandidateTierA == TRUE])
tierB <- order_candidates(x[CandidateTierB == TRUE])
tierC <- order_candidates(x[CandidateTierC == TRUE])

established <- order_candidates(
  x[EstablishedBenchmark == TRUE & CardiovascularATCGroup == TRUE]
)
# ------------------------------------------------------------
# Manuscript-facing compact table
# Tier A is intentionally small enough for manual literature review.
# ------------------------------------------------------------

compact_cols <- c(
  "Gene", "NatMed_DrugCategory", "NatMed_CategoryCode",
  "UniqueSNPs", "Tier1SNPs",
  "UniqueMetabolites", "Tier1Metabolites",
  "CodingSNPs", "SpliceRelatedSNPs",
  "HighImpactSNPs", "ModerateImpactSNPs",
  "OverlappingGeneSNPs", "NearestGeneSNPs",
  "MinMRFDR", "MaxPPH4", "MinM2P",
  "CPICTrackedPair", "CPICLevelABSupport",
  "DGIdbSupport",
  "ChEMBLMoASupport", "ChEMBLSingleProteinMoASupport",
  "ExternalDatabasePattern",
  "SNPs", "Metabolites",
  "CandidateClass"
)

compact_cols <- compact_cols[compact_cols %in% names(x)]

# ------------------------------------------------------------
# QC / selection flow
# ------------------------------------------------------------

qc <- data.table(
  Step = c(
    "All integrated gene-drug modules",
    "Cardiovascular drug-category modules",
    "CV modules with >=2 Tier-1 SNPs",
    "CV modules with >=2 Tier-1 SNPs and >=5 Tier-1 metabolites",
    "After excluding established CPIC A/B or ChEMBL single-protein modules",
    "Tier B: plus functional variant support",
    "Tier A: plus DGIdb support"
  ),
  N = c(
    nrow(x),
    x[CardiovascularATCGroup == TRUE, .N],
    x[CardiovascularATCGroup == TRUE & Tier1SNPs >= 2L, .N],
    x[InternalConvergence == TRUE, .N],
    x[CandidateTierC == TRUE, .N],
    x[CandidateTierB == TRUE, .N],
    x[CandidateTierA == TRUE, .N]
  )
)

class_qc <- x[, .N, by = CandidateClass][order(match(
  CandidateClass,
  c(
    "Established benchmark",
    "Tier A: convergent + functional + DGIdb",
    "Tier B: convergent + functional",
    "Tier C: convergent internal evidence",
    "Not prioritized"
  )
))]

# Drug-category distribution among Tier A candidates.
tierA_category_qc <- tierA[, .N, by = .(
  NatMed_CategoryCode, NatMed_DrugCategory
)][order(-N, NatMed_DrugCategory)]

# Gene-level recurrence across Tier A drug categories.
tierA_gene_qc <- tierA[, .(
  CandidateDrugCategories = uniqueN(NatMed_CategoryCode),
  DrugCategories = paste(sort(unique(NatMed_DrugCategory)), collapse = "; "),
  MaxTier1SNPs = max(Tier1SNPs, na.rm = TRUE),
  MaxTier1Metabolites = max(Tier1Metabolites, na.rm = TRUE),
  AnyCodingSupport = any(CodingSNPs > 0L, na.rm = TRUE),
  AnySpliceSupport = any(SpliceRelatedSNPs > 0L, na.rm = TRUE)
), by = Gene][order(-CandidateDrugCategories, -MaxTier1Metabolites, -MaxTier1SNPs, Gene)]

# ------------------------------------------------------------
# Write outputs
# ------------------------------------------------------------

fwrite(
  x,
  file.path(OUT_DIR, "CandidatePrioritization_AllModules.tsv"),
  sep = "\t", quote = FALSE, na = "NA"
)
fwrite(
  tierA,
  file.path(OUT_DIR, "CandidatePrioritization_TierA.tsv"),
  sep = "\t", quote = FALSE, na = "NA"
)
fwrite(
  tierB,
  file.path(OUT_DIR, "CandidatePrioritization_TierB.tsv"),
  sep = "\t", quote = FALSE, na = "NA"
)
fwrite(
  tierC,
  file.path(OUT_DIR, "CandidatePrioritization_TierC.tsv"),
  sep = "\t", quote = FALSE, na = "NA"
)
fwrite(
  established,
  file.path(OUT_DIR, "CandidatePrioritization_EstablishedBenchmarks.tsv"),
  sep = "\t", quote = FALSE, na = "NA"
)
fwrite(
  tierA[, ..compact_cols],
  file.path(OUT_DIR, "Table_CandidatePrioritization_TierA.tsv"),
  sep = "\t", quote = FALSE, na = "NA"
)
fwrite(
  qc,
  file.path(OUT_DIR, "CandidatePrioritization_QC.tsv"),
  sep = "\t", quote = FALSE, na = "NA"
)
fwrite(
  class_qc,
  file.path(OUT_DIR, "CandidatePrioritization_ClassQC.tsv"),
  sep = "\t", quote = FALSE, na = "NA"
)
fwrite(
  tierA_category_qc,
  file.path(OUT_DIR, "CandidatePrioritization_TierA_DrugCategoryQC.tsv"),
  sep = "\t", quote = FALSE, na = "NA"
)
fwrite(
  tierA_gene_qc,
  file.path(OUT_DIR, "CandidatePrioritization_TierA_GeneQC.tsv"),
  sep = "\t", quote = FALSE, na = "NA"
)

cat("\n============================================================\n")
cat("Candidate prioritization completed\n")
cat("============================================================\n\n")
print(qc)

cat("\nTier A candidates: ", nrow(tierA), "\n", sep = "")
cat("Tier B candidates: ", nrow(tierB), "\n", sep = "")
cat("Tier C candidates: ", nrow(tierC), "\n", sep = "")
cat("Established cardiovascular benchmark modules: ",
    nrow(established), "\n", sep = "")

cat("\nTier A preview:\n")
preview_cols <- c(
  "Gene", "NatMed_DrugCategory",
  "UniqueSNPs", "Tier1SNPs",
  "UniqueMetabolites", "Tier1Metabolites",
  "CodingSNPs", "SpliceRelatedSNPs",
  "ExternalDatabasePattern"
)
print(tierA[, ..preview_cols])

cat("\nOutputs:\n  ", OUT_DIR, "\n", sep = "")
cat("\nNEXT STEP:\n")
cat("Literature-review Tier A candidates before using terms such as\n")
cat("'less-established' or 'novel'. Absence from CPIC A/B or ChEMBL\n")
cat("single-protein MoA is not evidence of novelty.\n")
