# ============================================================
# Framework_v1
# Nature Medicine pharmacometabolomics functions
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
})

# ------------------------------------------------------------
# Text normalization for Nightingale trait matching
# ------------------------------------------------------------

normalize_trait_key <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x)] <- ""
  x <- iconv(x, from = "", to = "ASCII//TRANSLIT")
  x[is.na(x)] <- ""
  tolower(gsub("[^A-Za-z0-9]+", "", x))
}

sanitize_filename_token <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x
}

# ------------------------------------------------------------
# Curated aliases needed when abbreviation and explanation
# do not exactly match the Tambets/UKBB trait labels.
# ------------------------------------------------------------

natmed_manual_aliases <- function() {
  data.table(
    CAD_Key = normalize_trait_key(c(
      "VLDL_size",
      "Unsaturation",
      "LDL_size",
      "GlycA",
      "Total_FC",
      "Total_C",
      "Total_CE",
      "Phosphatidylc",
      "Phosphoglyc"
    )),
    NatMedAbbreviation = c(
      "VLDL.D",
      "UnsatDeg",
      "LDL.D",
      "Gp",
      "FreeC",
      "Serum.C",
      "EstC",
      "PC",
      "TotPG"
    )
  )
}

# ------------------------------------------------------------
# Read Suppl_T1: drug categories
# ------------------------------------------------------------

read_natmed_drug_categories <- function(workbook, sheet = "Suppl_T1") {
  if (!file.exists(workbook)) {
    stop_framework(paste0("Nature Medicine workbook not found: ", workbook))
  }
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop_framework("R package 'readxl' is required for the Nature Medicine workbook.")
  }

  x <- data.table::as.data.table(
    readxl::read_excel(
      workbook,
      sheet = sheet,
      skip = 1,
      col_types = "text"
    )
  )

  expected <- c(
    "Category code",
    "Drug name",
    "Anatomical main group",
    "Group letter",
    "Drug category",
    "ATC included"
  )

  missing <- setdiff(expected, names(x))
  if (length(missing) > 0) {
    stop_framework(
      paste0(
        "Suppl_T1 is missing expected columns: ",
        paste(missing, collapse = ", ")
      )
    )
  }

  x <- x[, .(
    CategoryCode = trimws(`Category code`),
    DrugName = trimws(`Drug name`),
    AnatomicalMainGroup = trimws(`Anatomical main group`),
    GroupLetter = trimws(`Group letter`),
    DrugCategory = trimws(`Drug category`),
    ATCIncluded = trimws(`ATC included`)
  )]

  x <- x[
    !is.na(DrugCategory) & DrugCategory != ""
  ]

  unique(x, by = "DrugCategory")
}

# ------------------------------------------------------------
# Read Suppl_T2: Nightingale metabolite dictionary
# ------------------------------------------------------------

read_natmed_metabolite_dictionary <- function(workbook, sheet = "Suppl_T2") {
  if (!file.exists(workbook)) {
    stop_framework(paste0("Nature Medicine workbook not found: ", workbook))
  }
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop_framework("R package 'readxl' is required for the Nature Medicine workbook.")
  }

  x <- data.table::as.data.table(
    readxl::read_excel(
      workbook,
      sheet = sheet,
      skip = 1,
      col_types = "text"
    )
  )

  expected <- c("Abbreviation", "Explanation", "Group")
  missing <- setdiff(expected, names(x))
  if (length(missing) > 0) {
    stop_framework(
      paste0(
        "Suppl_T2 is missing expected columns: ",
        paste(missing, collapse = ", ")
      )
    )
  }

  x <- x[, .(
    NatMedAbbreviation = trimws(Abbreviation),
    NatMedExplanation = trimws(Explanation),
    MetaboliteGroup = trimws(Group)
  )]

  x <- x[
    !is.na(NatMedAbbreviation) & NatMedAbbreviation != ""
  ]

  x[, AbbreviationKey := normalize_trait_key(NatMedAbbreviation)]
  x[, ExplanationKey := normalize_trait_key(NatMedExplanation)]

  unique(x, by = "NatMedAbbreviation")
}

# ------------------------------------------------------------
# Read Suppl_T4 M1/M2 association columns by position.
# Rows 1-3 are title/model/header rows, therefore skip = 3.
# ------------------------------------------------------------

read_natmed_m1_m2_associations <- function(workbook, sheet = "Suppl_T4") {
  if (!file.exists(workbook)) {
    stop_framework(paste0("Nature Medicine workbook not found: ", workbook))
  }
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop_framework("R package 'readxl' is required for the Nature Medicine workbook.")
  }

  raw <- data.table::as.data.table(
    readxl::read_excel(
      workbook,
      sheet = sheet,
      skip = 3,
      col_names = FALSE,
      col_types = "text"
    )
  )

  if (ncol(raw) < 16) {
    stop_framework(
      paste0(
        "Suppl_T4 contains only ", ncol(raw),
        " columns after import; at least 16 were expected."
      )
    )
  }

  x <- raw[, 1:16]
  setnames(
    x,
    c(
      "NatMedAbbreviation",
      "DrugCategory",
      "M1_UserN",
      "M1_NonUserN",
      "M1_Beta",
      "M1_SE",
      "M1_P",
      "M1_Threshold",
      "M1_SignificanceFlag",
      "M2_UserN",
      "M2_NonUserN",
      "M2_Beta",
      "M2_SE",
      "M2_P",
      "M2_Threshold",
      "M1M2_SignificanceFlag"
    )
  )

  x[, NatMedAbbreviation := trimws(as.character(NatMedAbbreviation))]
  x[, DrugCategory := trimws(as.character(DrugCategory))]

  numeric_columns <- c(
    "M1_UserN", "M1_NonUserN", "M1_Beta", "M1_SE", "M1_P", "M1_Threshold",
    "M2_UserN", "M2_NonUserN", "M2_Beta", "M2_SE", "M2_P", "M2_Threshold"
  )

  for (column in numeric_columns) {
    set(
      x,
      j = column,
      value = suppressWarnings(as.numeric(x[[column]]))
    )
  }

  x <- x[
    !is.na(NatMedAbbreviation) & NatMedAbbreviation != "" &
      !is.na(DrugCategory) & DrugCategory != ""
  ]

  x[, M1Significant :=
      !is.na(M1_P) & !is.na(M1_Threshold) & M1_P < M1_Threshold]

  x[, M2Significant :=
      !is.na(M2_P) & !is.na(M2_Threshold) & M2_P < M2_Threshold]

  # The Suppl_T4 column explicitly means significant in both M1 and M2.
  x[, M1M2Significant := trimws(M1M2_SignificanceFlag) == "*"]

  x[, M1M2DirectionConsistent :=
      !is.na(M1_Beta) & !is.na(M2_Beta) &
        sign(M1_Beta) == sign(M2_Beta)]

  x[, M1M2Robust := M1M2Significant & M1M2DirectionConsistent]

  x
}

# ------------------------------------------------------------
# Harmonize CAD/Tambets metabolite traits to NatMed abbreviations
# ------------------------------------------------------------

map_coloc_to_natmed <- function(coloc, metabolite_dictionary) {
  coloc <- data.table::copy(coloc)
  dict <- data.table::copy(metabolite_dictionary)

  required_coloc <- c(
    "metabolite",
    "study_accession",
    "PP.H4",
    "reported_trait",
    "ftp_path"
  )
  check_columns(coloc, required_coloc)

  coloc[, CAD_MetaboliteKey := normalize_trait_key(metabolite)]
  coloc[, CAD_ReportedTraitKey := normalize_trait_key(reported_trait)]

  dict_abbr <- dict[, .(
    AbbreviationKey,
    NatMedAbbreviation_abbr = NatMedAbbreviation,
    NatMedExplanation_abbr = NatMedExplanation,
    MetaboliteGroup_abbr = MetaboliteGroup
  )]

  result <- merge(
    coloc,
    dict_abbr,
    by.x = "CAD_MetaboliteKey",
    by.y = "AbbreviationKey",
    all.x = TRUE,
    sort = FALSE
  )

  result[, RowID := .I]

  result[, `:=`(
    NatMedAbbreviation = NatMedAbbreviation_abbr,
    NatMedExplanation = NatMedExplanation_abbr,
    MetaboliteGroup = MetaboliteGroup_abbr,
    MatchMethod = fifelse(!is.na(NatMedAbbreviation_abbr), "abbreviation", NA_character_)
  )]

  result[, c(
    "NatMedAbbreviation_abbr",
    "NatMedExplanation_abbr",
    "MetaboliteGroup_abbr"
  ) := NULL]

  # Match unresolved records by exact normalized reported-trait explanation.
  unresolved <- is.na(result$NatMedAbbreviation)
  if (any(unresolved)) {
    dict_exp <- dict[, .(
      ExplanationKey,
      NatMedAbbreviation_exp = NatMedAbbreviation,
      NatMedExplanation_exp = NatMedExplanation,
      MetaboliteGroup_exp = MetaboliteGroup
    )]

    exp_matches <- merge(
      result[unresolved, .(RowID, CAD_ReportedTraitKey)],
      dict_exp,
      by.x = "CAD_ReportedTraitKey",
      by.y = "ExplanationKey",
      all.x = TRUE,
      sort = FALSE
    )

    hit <- !is.na(exp_matches$NatMedAbbreviation_exp)
    if (any(hit)) {
      target <- exp_matches$RowID[hit]
      result[target, `:=`(
        NatMedAbbreviation = exp_matches$NatMedAbbreviation_exp[hit],
        NatMedExplanation = exp_matches$NatMedExplanation_exp[hit],
        MetaboliteGroup = exp_matches$MetaboliteGroup_exp[hit],
        MatchMethod = "reported_trait_explanation"
      )]
    }
  }

  # Apply a small, explicit alias table only to remaining unresolved traits.
  aliases <- natmed_manual_aliases()
  alias_dict <- merge(
    aliases,
    dict[, .(
      NatMedAbbreviation,
      NatMedExplanation,
      MetaboliteGroup
    )],
    by = "NatMedAbbreviation",
    all.x = TRUE,
    sort = FALSE
  )

  unresolved <- is.na(result$NatMedAbbreviation)
  if (any(unresolved)) {
    alias_matches <- merge(
      result[unresolved, .(RowID, CAD_MetaboliteKey)],
      alias_dict,
      by.x = "CAD_MetaboliteKey",
      by.y = "CAD_Key",
      all.x = TRUE,
      sort = FALSE
    )

    hit <- !is.na(alias_matches$NatMedAbbreviation)
    if (any(hit)) {
      target <- alias_matches$RowID[hit]
      result[target, `:=`(
        NatMedAbbreviation = alias_matches$NatMedAbbreviation[hit],
        NatMedExplanation = alias_matches$NatMedExplanation[hit],
        MetaboliteGroup = alias_matches$MetaboliteGroup[hit],
        MatchMethod = "curated_alias"
      )]
    }
  }

  result[, MatchStatus := fifelse(
    !is.na(NatMedAbbreviation) & NatMedAbbreviation != "",
    "matched",
    "not_in_NatMed_150_trait_panel"
  )]

  result[, RowID := NULL]
  result[]
}

# ------------------------------------------------------------
# Build pharmacometabolomic evidence by joining MR with M2
# ------------------------------------------------------------

classify_drug_metabolite_direction <- function(beta) {
  fifelse(
    is.na(beta),
    NA_character_,
    fifelse(beta > 0, "Higher_in_drug_users",
      fifelse(beta < 0, "Lower_in_drug_users", "Neutral")
    )
  )
}

calculate_directional_concordance <- function(mr_beta, drug_beta) {
  valid <- !is.na(mr_beta) & !is.na(drug_beta) & mr_beta != 0 & drug_beta != 0
  out <- rep(NA, length(mr_beta))
  out[valid] <- sign(mr_beta[valid]) != sign(drug_beta[valid])
  out
}

build_pharmacometabolomic_evidence <- function(
    mr_summary,
    candidate_mapping,
    m2_associations,
    drug_categories,
    fallback_threshold = 1.9157088122605367e-5
) {
  mr <- data.table::copy(mr_summary)
  mapping <- data.table::copy(candidate_mapping)
  assoc <- data.table::copy(m2_associations)
  drugs <- data.table::copy(drug_categories)

  check_columns(
    mr,
    c(
      "Metabolite", "MR_beta", "MR_se", "MR_p", "MR_FDR",
      "PP_H4", "OR", "OR_lci95", "OR_uci95", "nsnp"
    )
  )

  check_columns(
    mapping,
    c(
      "metabolite", "NatMedAbbreviation", "MatchStatus",
      "NatMedExplanation", "MetaboliteGroup",
      "reported_trait", "study_accession", "PP.H4", "MatchMethod"
    )
  )

  check_columns(
    assoc,
    c(
      "NatMedAbbreviation", "DrugCategory",
      "M1_Beta", "M1_P", "M2_Beta", "M2_P"
    )
  )

  check_columns(
    drugs,
    c(
      "DrugCategory", "DrugName", "CategoryCode",
      "AnatomicalMainGroup", "GroupLetter", "ATCIncluded"
    )
  )

  # ----------------------------------------------------------
  # One canonical mapping row per CAD/Nightingale trait.
  # ----------------------------------------------------------

  map_use <- unique(
    mapping[MatchStatus == "matched", .(
      Metabolite = metabolite,
      NatMedAbbreviation,
      NatMedExplanation,
      MetaboliteGroup,
      reported_trait,
      study_accession,
      PP_H4_coloc = as.numeric(PP.H4),
      MatchMethod
    )],
    by = "Metabolite"
  )

  if (anyDuplicated(map_use$Metabolite)) {
    stop_framework("Duplicate metabolite rows remain in the NatMed mapping.")
  }

  if (anyDuplicated(map_use$NatMedAbbreviation)) {
    stop_framework(
      paste0(
        "More than one CAD trait maps to the same NatMed abbreviation. ",
        "This must be reviewed before pharmacometabolomic integration."
      )
    )
  }

  # Keep only the 102 MR traits that passed NatMed mapping.
  mr <- merge(
    mr,
    map_use,
    by = "Metabolite",
    all = FALSE,
    sort = FALSE
  )

  # ----------------------------------------------------------
  # M2 significance: use the row-specific published threshold
  # when present; otherwise use the published fallback.
  # ----------------------------------------------------------

  if (!"M2_Threshold" %in% names(assoc)) {
    assoc[, M2_Threshold := NA_real_]
  }

  assoc[, EffectiveM2Threshold := fifelse(
    !is.na(M2_Threshold),
    as.numeric(M2_Threshold),
    fallback_threshold
  )]

  assoc[, M2Significant :=
    !is.na(M2_P) &
      !is.na(EffectiveM2Threshold) &
      as.numeric(M2_P) < EffectiveM2Threshold
  ]

  # ----------------------------------------------------------
  # IMPORTANT:
  # 00_prepare_natmed_inputs.R already enriches Suppl_T4 with
  # metabolite dictionary and drug metadata. To avoid duplicate
  # .x/.y columns, retain only association/model columns here
  # and add canonical metadata exactly once below.
  # ----------------------------------------------------------

  association_columns <- intersect(
    c(
      "NatMedAbbreviation",
      "DrugCategory",
      "M1_UserN",
      "M1_NonUserN",
      "M1_Beta",
      "M1_SE",
      "M1_P",
      "M1_Threshold",
      "M1_SignificanceFlag",
      "M1Significant",
      "M2_UserN",
      "M2_NonUserN",
      "M2_Beta",
      "M2_SE",
      "M2_P",
      "M2_Threshold",
      "M1M2_SignificanceFlag",
      "M1M2Significant",
      "M1M2DirectionConsistent",
      "M1M2Robust",
      "EffectiveM2Threshold",
      "M2Significant"
    ),
    names(assoc)
  )

  assoc_core <- assoc[, ..association_columns]

  # Guard against unexpected duplicate metabolite-drug rows.
  duplicate_edges <- assoc_core[
    ,
    .N,
    by = .(NatMedAbbreviation, DrugCategory)
  ][N > 1]

  if (nrow(duplicate_edges) > 0) {
    stop_framework(
      paste0(
        "Duplicate NatMed metabolite-drug rows detected in Suppl_T4: ",
        nrow(duplicate_edges)
      )
    )
  }

  # Expand each MR metabolite to its 87 NatMed drug-category rows.
  result <- merge(
    mr,
    assoc_core,
    by = "NatMedAbbreviation",
    allow.cartesian = TRUE,
    all.x = TRUE,
    sort = FALSE
  )

  # Add drug metadata exactly once.
  drug_meta <- unique(
    drugs[, .(
      DrugCategory,
      CategoryCode,
      DrugName,
      AnatomicalMainGroup,
      GroupLetter,
      ATCIncluded
    )],
    by = "DrugCategory"
  )

  result <- merge(
    result,
    drug_meta,
    by = "DrugCategory",
    all.x = TRUE,
    sort = FALSE
  )

  if (any(is.na(result$CategoryCode))) {
    stop_framework(
      "One or more pharmacometabolomic rows failed to map to Suppl_T1 drug metadata."
    )
  }

  # ----------------------------------------------------------
  # Directional interpretation
  #
  # MR_beta > 0: genetically higher metabolite -> higher CAD risk
  # M2_Beta < 0: drug users have lower metabolite level
  #
  # Opposite beta signs are therefore labelled risk-opposing.
  # This is NOT interpreted as a causal drug effect.
  # ----------------------------------------------------------

  result[, MR_Direction := fifelse(
    MR_beta > 0,
    "Higher_metabolite_higher_CAD_risk",
    fifelse(
      MR_beta < 0,
      "Higher_metabolite_lower_CAD_risk",
      "Neutral"
    )
  )]

  result[, DrugMetaboliteDirection :=
    classify_drug_metabolite_direction(M2_Beta)
  ]

  result[, RiskOpposingDirection :=
    calculate_directional_concordance(MR_beta, M2_Beta)
  ]

  # Backward-compatible field name; explicitly means
  # therapeutic/risk-opposing concordance, not equal beta signs.
  result[, DirectionalConcordance := RiskOpposingDirection]

  result[, DirectionClass := fifelse(
    is.na(RiskOpposingDirection),
    "Indeterminate",
    fifelse(
      RiskOpposingDirection,
      "Risk_opposing",
      "Risk_reinforcing"
    )
  )]

  result[, CardiovascularATCGroup :=
    !is.na(GroupLetter) & GroupLetter == "C"
  ]

  result[, MR_FDRSignificant :=
    !is.na(MR_FDR) & MR_FDR < 0.05
  ]

  result[, PrimaryPharmacoGeneticEvidence :=
    M2Significant %in% TRUE &
      RiskOpposingDirection %in% TRUE &
      MR_FDRSignificant %in% TRUE
  ]

  result[, EvidenceTier := fifelse(
    PrimaryPharmacoGeneticEvidence %in% TRUE &
      M1M2Robust %in% TRUE,
    "Tier1_Robust",
    fifelse(
      PrimaryPharmacoGeneticEvidence %in% TRUE,
      "Tier2_HighConfidence",
      fifelse(
        M2Significant %in% TRUE &
          MR_FDRSignificant %in% TRUE &
          RiskOpposingDirection %in% FALSE,
        "DirectionallyDiscordant",
        "NotPrimary"
      )
    )
  )]

  result[, PharmacoEvidenceType :=
    "Cross-sectional drug-metabolite association (NatMed M2)"
  ]

  preferred <- c(
    "Metabolite",
    "NatMedAbbreviation",
    "NatMedExplanation",
    "MetaboliteGroup",
    "reported_trait",
    "study_accession",
    "PP_H4",
    "PP_H4_coloc",
    "MatchMethod",
    "MR_beta",
    "MR_se",
    "MR_p",
    "MR_FDR",
    "MR_FDRSignificant",
    "OR",
    "OR_lci95",
    "OR_uci95",
    "nsnp",
    "MR_Direction",
    "DrugCategory",
    "CategoryCode",
    "DrugName",
    "AnatomicalMainGroup",
    "GroupLetter",
    "ATCIncluded",
    "CardiovascularATCGroup",
    "M1_UserN",
    "M1_NonUserN",
    "M1_Beta",
    "M1_SE",
    "M1_P",
    "M1_Threshold",
    "M1Significant",
    "M2_UserN",
    "M2_NonUserN",
    "M2_Beta",
    "M2_SE",
    "M2_P",
    "M2_Threshold",
    "EffectiveM2Threshold",
    "M2Significant",
    "M1M2Significant",
    "M1M2DirectionConsistent",
    "M1M2Robust",
    "DrugMetaboliteDirection",
    "RiskOpposingDirection",
    "DirectionalConcordance",
    "DirectionClass",
    "PrimaryPharmacoGeneticEvidence",
    "EvidenceTier",
    "PharmacoEvidenceType"
  )

  keep_cols <- intersect(preferred, names(result))
  result <- result[, ..keep_cols]
  setorder(result, Metabolite, M2_P, DrugCategory)

  result[]
}
