# ============================================================
# Framework_v1
# Evidence-matrix functions -- Nature Medicine version
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
})

classify_mr_evidence <- function(mr_fdr, mr_p,
                                 fdr_threshold = 0.05,
                                 nominal_threshold = 0.05) {
  fcase(
    !is.na(mr_fdr) & mr_fdr < fdr_threshold, "FDRSignificant",
    !is.na(mr_p) & mr_p < nominal_threshold, "NominalEvidence",
    default = "NoMREvidence"
  )
}

classify_evidence_tier <- function(
    m2_significant,
    directional_concordance,
    mr_fdr,
    mr_p,
    m1m2_robust,
    fdr_threshold = 0.05,
    nominal_threshold = 0.05
) {
  fcase(
    m2_significant %in% TRUE &
      directional_concordance %in% TRUE &
      !is.na(mr_fdr) & mr_fdr < fdr_threshold &
      m1m2_robust %in% TRUE,
    "Tier1_Robust",

    m2_significant %in% TRUE &
      directional_concordance %in% TRUE &
      !is.na(mr_fdr) & mr_fdr < fdr_threshold,
    "Tier2_HighConfidence",

    m2_significant %in% TRUE &
      directional_concordance %in% TRUE &
      !is.na(mr_p) & mr_p < nominal_threshold,
    "Tier3_NominalMR",

    m2_significant %in% TRUE &
      directional_concordance %in% TRUE,
    "Tier4_Pharmacometabolic",

    m2_significant %in% TRUE &
      directional_concordance %in% FALSE,
    "Discordant",

    m2_significant %in% TRUE,
    "DirectionIndeterminate",

    default = "M2_NotSignificant"
  )
}

build_evidence_matrix <- function(
    pharmacometabolomics,
    mr_fdr_threshold = 0.05,
    mr_p_nominal_threshold = 0.05,
    require_directional_concordance = TRUE,
    require_mr_fdr = TRUE
) {
  x <- data.table::copy(pharmacometabolomics)

  check_columns(
    x,
    c(
      "Metabolite",
      "DrugCategory",
      "MR_beta",
      "MR_p",
      "MR_FDR",
      "PP_H4",
      "M2_Beta",
      "M2_P",
      "M2Significant",
      "DirectionalConcordance",
      "M1M2Robust"
    )
  )

  x[, MREvidenceClass := classify_mr_evidence(
    MR_FDR,
    MR_p,
    fdr_threshold = mr_fdr_threshold,
    nominal_threshold = mr_p_nominal_threshold
  )]

  x[, ColocalizationEvidenceClass := fifelse(
    !is.na(PP_H4) & PP_H4 >= 0.99,
    "VeryStrong_PP.H4>=0.99",
    fifelse(
      !is.na(PP_H4),
      "BelowPrimaryThreshold",
      "Unavailable"
    )
  )]

  x[, PharmacometabolomicEvidenceClass := fifelse(
    M2Significant,
    "M2_MultipleTestingSignificant",
    "M2_NotSignificant"
  )]

  x[, EvidenceTier := classify_evidence_tier(
    m2_significant = M2Significant,
    directional_concordance = DirectionalConcordance,
    mr_fdr = MR_FDR,
    mr_p = MR_p,
    m1m2_robust = M1M2Robust,
    fdr_threshold = mr_fdr_threshold,
    nominal_threshold = mr_p_nominal_threshold
  )]

  primary <- x$M2Significant %in% TRUE

  if (isTRUE(require_directional_concordance)) {
    primary <- primary & x$DirectionalConcordance %in% TRUE
  }

  if (isTRUE(require_mr_fdr)) {
    primary <- primary & !is.na(x$MR_FDR) & x$MR_FDR < mr_fdr_threshold
  }

  x[, CandidateSelected := primary]

  x[, CandidateSelectionReason := fifelse(
    CandidateSelected & M1M2Robust,
    "M2Significant;OpposesMRRiskDirection;MRFDRSignificant;M1M2Robust",
    fifelse(
      CandidateSelected,
      "M2Significant;OpposesMRRiskDirection;MRFDRSignificant",
      NA_character_
    )
  )]

  # No arbitrary composite score is used. Ranking is lexicographic and
  # transparently based on evidence tier, M2 P value, MR FDR and effect size.
  tier_order <- c(
    "Tier1_Robust" = 1L,
    "Tier2_HighConfidence" = 2L,
    "Tier3_NominalMR" = 3L,
    "Tier4_Pharmacometabolic" = 4L,
    "DirectionIndeterminate" = 5L,
    "Discordant" = 6L,
    "M2_NotSignificant" = 7L
  )

  x[, TierOrder := unname(tier_order[EvidenceTier])]
  x[is.na(TierOrder), TierOrder := 99L]

  x[, M2PForRank := fifelse(is.na(M2_P), Inf, M2_P)]
  x[, MRFDRForRank := fifelse(is.na(MR_FDR), Inf, MR_FDR)]
  x[, AbsM2BetaForRank := fifelse(is.na(M2_Beta), -Inf, abs(M2_Beta))]

  x[, MetabolitePriorityRank := frankv(
    .SD,
    cols = c("TierOrder", "M2PForRank", "MRFDRForRank", "AbsM2BetaForRank"),
    order = c(1L, 1L, 1L, -1L),
    ties.method = "min"
  ), by = Metabolite,
  .SDcols = c("TierOrder", "M2PForRank", "MRFDRForRank", "AbsM2BetaForRank")]

  x[, GlobalPriorityRank := frankv(
    x,
    cols = c("TierOrder", "M2PForRank", "MRFDRForRank", "AbsM2BetaForRank"),
    order = c(1L, 1L, 1L, -1L),
    ties.method = "min"
  )]

  x[, c("M2PForRank", "MRFDRForRank", "AbsM2BetaForRank") := NULL]

  data.table::setorderv(
    x,
    cols = c("CandidateSelected", "TierOrder", "MetabolitePriorityRank", "M2_P", "DrugCategory"),
    order = c(-1L, 1L, 1L, 1L, 1L),
    na.last = TRUE
  )

  x
}

select_drug_class_candidates <- function(evidence_matrix) {
  candidates <- copy(evidence_matrix)[CandidateSelected %in% TRUE]
  setorder(candidates, GlobalPriorityRank, Metabolite, MetabolitePriorityRank)
  candidates
}

summarize_evidence_matrix <- function(evidence_matrix, candidate_network) {
  data.table(
    Metric = c(
      "Evidence matrix rows",
      "Metabolites in evidence matrix",
      "Drug categories in evidence matrix",
      "M2-significant associations",
      "Directionally concordant M2-significant associations",
      "Primary candidate associations",
      "Primary candidate metabolites",
      "Primary candidate drug categories",
      "Tier 1 primary candidates",
      "Tier 2 primary candidates"
    ),
    Value = c(
      nrow(evidence_matrix),
      uniqueN(evidence_matrix$Metabolite),
      uniqueN(evidence_matrix$DrugCategory),
      evidence_matrix[M2Significant %in% TRUE, .N],
      evidence_matrix[M2Significant %in% TRUE & DirectionalConcordance %in% TRUE, .N],
      nrow(candidate_network),
      uniqueN(candidate_network$Metabolite),
      uniqueN(candidate_network$DrugCategory),
      candidate_network[EvidenceTier == "Tier1_Robust", .N],
      candidate_network[EvidenceTier == "Tier2_HighConfidence", .N]
    )
  )
}

summarize_candidates_by_metabolite <- function(candidate_network) {
  candidate_network[
    ,
    .(
      N_CandidateDrugClasses = .N,
      N_Tier1 = sum(EvidenceTier == "Tier1_Robust", na.rm = TRUE),
      N_Tier2 = sum(EvidenceTier == "Tier2_HighConfidence", na.rm = TRUE),
      BestM2P = suppressWarnings(min(M2_P, na.rm = TRUE)),
      BestDrugCategory = DrugCategory[which.min(M2_P)][1]
    ),
    by = .(
      Metabolite,
      NatMedAbbreviation,
      NatMedExplanation,
      MetaboliteGroup,
      MR_beta,
      MR_p,
      MR_FDR,
      PP_H4
    )
  ]
}
