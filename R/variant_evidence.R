# ============================================================
# Framework_v1
# Variant-level evidence functions
# ============================================================


# ------------------------------------------------------------
# Resolve a column from possible aliases
# ------------------------------------------------------------

resolve_column_alias <- function(
    data,
    aliases,
    required = FALSE,
    label = NULL
){
  
  matched <- aliases[
    aliases %in% names(data)
  ]
  
  if(length(matched) > 0){
    
    return(
      matched[1]
    )
  }
  
  if(isTRUE(required)){
    
    if(is.null(label)){
      
      label <- paste(
        aliases,
        collapse = " / "
      )
    }
    
    stop_framework(
      paste0(
        "Required column could not be resolved: ",
        label,
        ". Available columns: ",
        paste(
          names(data),
          collapse = ", "
        )
      )
    )
  }
  
  NA_character_
}


# ------------------------------------------------------------
# Copy an aliased column to a standard column name
# ------------------------------------------------------------

standardize_alias_column <- function(
    data,
    standard_name,
    aliases,
    required = FALSE
){
  
  source_column <- resolve_column_alias(
    data = data,
    aliases = unique(
      c(
        standard_name,
        aliases
      )
    ),
    required = required,
    label = standard_name
  )
  
  if(is.na(source_column)){
    
    if(!standard_name %in% names(data)){
      
      data[
        ,
        (standard_name) := NA
      ]
    }
    
    return(
      invisible(data)
    )
  }
  
  if(source_column != standard_name){
    
    data[
      ,
      (standard_name) :=
        get(source_column)
    ]
  }
  
  invisible(data)
}


# ------------------------------------------------------------
# Convert likely logical fields
# ------------------------------------------------------------

as_framework_logical <- function(x){
  
  if(is.logical(x)){
    
    return(x)
  }
  
  normalized <- toupper(
    trimws(
      as.character(x)
    )
  )
  
  result <- rep(
    NA,
    length(normalized)
  )
  
  result[
    normalized %in% c(
      "TRUE",
      "T",
      "1",
      "YES",
      "Y"
    )
  ] <- TRUE
  
  result[
    normalized %in% c(
      "FALSE",
      "F",
      "0",
      "NO",
      "N"
    )
  ] <- FALSE
  
  result
}


# ------------------------------------------------------------
# Read and standardize harmonized instruments
# ------------------------------------------------------------

read_harmonized_instruments <- function(file){
  
  if(!file.exists(file)){
    
    stop_framework(
      paste0(
        "Harmonized instrument file not found: ",
        file
      )
    )
  }
  
  instruments <- data.table::as.data.table(
    read_tsv_safe(file)
  )
  
  # --------------------------------------------------------
  # Required identifiers
  # --------------------------------------------------------
  
  standardize_alias_column(
    instruments,
    standard_name = "Metabolite",
    aliases = c(
      "metabolite",
      "Exposure",
      "exposure",
      "exposure_name",
      "trait",
      "Trait"
    ),
    required = FALSE
  )
  
  standardize_alias_column(
    instruments,
    standard_name = "SNP",
    aliases = c(
      "rsid",
      "RSID",
      "MarkerName",
      "variant",
      "Variant",
      "variant_id",
      "ID"
    ),
    required = TRUE
  )
  
  # --------------------------------------------------------
  # Study identifiers
  # --------------------------------------------------------
  
  standardize_alias_column(
    instruments,
    standard_name = "StudyAccession",
    aliases = c(
      "study_accession",
      "Study",
      "study",
      "ExposureStudy",
      "id.exposure",
      "exposure_id"
    )
  )
  
  standardize_alias_column(
    instruments,
    standard_name = "OutcomeStudy",
    aliases = c(
      "outcome",
      "Outcome",
      "outcome_name"
    )
  )
  # --------------------------------------------------------
  # Genomic position
  # --------------------------------------------------------
  
  standardize_alias_column(
    instruments,
    standard_name = "CHR",
    aliases = c(
      "chr",
      "chromosome",
      "Chromosome",
      "CHR_exposure",
      "chr.exposure"
    )
  )
  
  standardize_alias_column(
    instruments,
    standard_name = "BP",
    aliases = c(
      "bp",
      "position",
      "Position",
      "POS",
      "pos",
      "BP_exposure",
      "pos.exposure"
    )
  )
  
  # --------------------------------------------------------
  # Alleles
  # --------------------------------------------------------
  
  standardize_alias_column(
    instruments,
    standard_name = "EffectAllele",
    aliases = c(
      "effect_allele",
      "EA",
      "A1",
      "Allele1",
      "effect_allele.exposure"
    )
  )
  
  standardize_alias_column(
    instruments,
    standard_name = "OtherAllele",
    aliases = c(
      "other_allele",
      "OA",
      "A2",
      "Allele2",
      "other_allele.exposure"
    )
  )
  
  standardize_alias_column(
    instruments,
    standard_name = "EffectAlleleFrequency",
    aliases = c(
      "eaf",
      "EAF",
      "effect_allele_frequency",
      "eaf.exposure"
    )
  )
  
  # --------------------------------------------------------
  # Exposure statistics
  # --------------------------------------------------------
  
  standardize_alias_column(
    instruments,
    standard_name = "BetaExposure",
    aliases = c(
      "beta_exposure",
      "Beta_exposure",
      "beta.exposure",
      "beta",
      "Effect"
    )
  )
  
  standardize_alias_column(
    instruments,
    standard_name = "SEExposure",
    aliases = c(
      "se_exposure",
      "SE_exposure",
      "se.exposure",
      "se",
      "StdErr"
    )
  )
  
  standardize_alias_column(
    instruments,
    standard_name = "PExposure",
    aliases = c(
      "p_exposure",
      "P_exposure",
      "pval.exposure",
      "pval",
      "P",
      "P_value",
      "P-value"
    )
  )
  
  # --------------------------------------------------------
  # Outcome statistics
  # --------------------------------------------------------
  
  standardize_alias_column(
    instruments,
    standard_name = "BetaOutcome",
    aliases = c(
      "beta_outcome",
      "Beta_outcome",
      "beta.outcome"
    )
  )
  
  standardize_alias_column(
    instruments,
    standard_name = "SEOutcome",
    aliases = c(
      "se_outcome",
      "SE_outcome",
      "se.outcome"
    )
  )
  
  standardize_alias_column(
    instruments,
    standard_name = "POutcome",
    aliases = c(
      "p_outcome",
      "P_outcome",
      "pval.outcome"
    )
  )
  
  standardize_alias_column(
    instruments,
    standard_name = "EffectAlleleOutcome",
    aliases = c(
      "effect_allele.outcome",
      "effect_allele_outcome"
    )
  )
  
  standardize_alias_column(
    instruments,
    standard_name = "OtherAlleleOutcome",
    aliases = c(
      "other_allele.outcome",
      "other_allele_outcome"
    )
  )
  
  standardize_alias_column(
    instruments,
    standard_name = "EffectAlleleFrequencyOutcome",
    aliases = c(
      "eaf.outcome",
      "eaf_outcome"
    )
  )
  
  standardize_alias_column(
    instruments,
    standard_name = "CHROutcome",
    aliases = c(
      "CHR_outcome",
      "chr.outcome",
      "chr_outcome"
    )
  )
  
  standardize_alias_column(
    instruments,
    standard_name = "BPOutcome",
    aliases = c(
      "BP_outcome",
      "pos.outcome",
      "bp_outcome"
    )
  )
  
  
  # --------------------------------------------------------
  # Instrument quality
  # --------------------------------------------------------
  
  standardize_alias_column(
    instruments,
    standard_name = "FStatistic",
    aliases = c(
      "F_statistic",
      "F",
      "Fstat",
      "F_stat",
      "instrument_F"
    )
  )
  
  standardize_alias_column(
    instruments,
    standard_name = "R2",
    aliases = c(
      "r2",
      "R_squared",
      "R2_exposure",
      "instrument_R2"
    )
  )
  
  standardize_alias_column(
    instruments,
    standard_name = "Palindromic",
    aliases = c(
      "palindromic",
      "PalindromicSNP"
    )
  )
  
  standardize_alias_column(
    instruments,
    standard_name = "Ambiguous",
    aliases = c(
      "ambiguous",
      "AmbiguousSNP"
    )
  )
  
  standardize_alias_column(
    instruments,
    standard_name = "HarmonizationStatus",
    aliases = c(
      "harmonisation",
      "harmonization",
      "harmonization_status",
      "Harmonised",
      "harmonized",
      "mr_keep",
      "MRKeep"
    )
  )
  
  standardize_alias_column(
    instruments,
    standard_name = "MHC",
    aliases = c(
      "MHC_status",
      "InMHC",
      "in_mhc",
      "MHC_flag"
    )
  )
  
  standardize_alias_column(
    instruments,
    standard_name = "SteigerDirection",
    aliases = c(
      "steiger_dir",
      "Steiger_status",
      "correct_causal_direction",
      "steiger_direction"
    )
  )
  
  # --------------------------------------------------------
  # Numeric conversion
  # --------------------------------------------------------
  
  numeric_columns <- intersect(
    c(
      "CHR",
      "BP",
      "CHROutcome",
      "BPOutcome",
      "EffectAlleleFrequency",
      "EffectAlleleFrequencyOutcome",
      "BetaExposure",
      "SEExposure",
      "PExposure",
      "BetaOutcome",
      "SEOutcome",
      "POutcome",
      "FStatistic",
      "R2"
    ),
    names(instruments)
  )
  
  for(column in numeric_columns){
    
    data.table::set(
      instruments,
      j = column,
      value = suppressWarnings(
        as.numeric(
          instruments[[column]]
        )
      )
    )
  }
  
  logical_columns <- intersect(
    c(
      "Palindromic",
      "Ambiguous",
      "MHC"
    ),
    names(instruments)
  )
  
  for(column in logical_columns){
    
    data.table::set(
      instruments,
      j = column,
      value = as_framework_logical(
        instruments[[column]]
      )
    )
  }
  
  instruments[
    ,
    Metabolite :=
      trimws(
        as.character(
          Metabolite
        )
      )
  ]
  
  instruments[
    ,
    SNP :=
      trimws(
        as.character(
          SNP
        )
      )
  ]
  
  instruments <- instruments[
    !is.na(SNP) &
      SNP != ""
  ]
  
  instruments
}


# ------------------------------------------------------------
# Read all metabolite-specific harmonized instrument files
# ------------------------------------------------------------

read_harmonized_instrument_directory <- function(
    directory,
    pattern = "_harmonized_for_MR_rsID\\.tsv$"
){
  
  if(!dir.exists(directory)){
    
    stop_framework(
      paste0(
        "Harmonized instrument directory not found: ",
        directory
      )
    )
  }
  
  files <- list.files(
    path = directory,
    pattern = pattern,
    full.names = TRUE
  )
  
  if(length(files) == 0){
    
    stop_framework(
      paste0(
        "No harmonized instrument files found in: ",
        directory,
        " using pattern: ",
        pattern
      )
    )
  }
  
  message(
    "Reading ",
    length(files),
    " harmonized instrument files."
  )
  
  tables <- lapply(
    files,
    function(file){
      
      metabolite_from_file <- sub(
        pattern,
        "",
        basename(file)
      )
      
      table <- read_harmonized_instruments(
        file
      )
      
      # Dosyada metabolit kolonu yoksa veya yanlışsa
      # dosya adındaki metabolit adını kullan.
      if(
        !"Metabolite" %in% names(table) ||
        all(
          is.na(table$Metabolite) |
          table$Metabolite == ""
        )
      ){
        
        table[
          ,
          Metabolite :=
            metabolite_from_file
        ]
      }
      
      table[
        ,
        SourceFile :=
          basename(file)
      ]
      
      table
    }
  )
  
  combined <- data.table::rbindlist(
    tables,
    use.names = TRUE,
    fill = TRUE
  )
  
  combined[
    ,
    Metabolite :=
      trimws(
        as.character(
          Metabolite
        )
      )
  ]
  
  combined
}


# ------------------------------------------------------------
# Complete missing F statistics where possible
# ------------------------------------------------------------

calculate_missing_instrument_strength <- function(
    instruments
){
  
  result <- data.table::copy(
    instruments
  )
  
  if(!"FStatistic" %in% names(result)){
    
    result[
      ,
      FStatistic := NA_real_
    ]
  }
  
  computable <- (
    is.na(result$FStatistic) |
      !is.finite(result$FStatistic)
  ) &
    !is.na(result$BetaExposure) &
    !is.na(result$SEExposure) &
    result$SEExposure > 0
  
  result[
    computable,
    FStatistic :=
      (
        BetaExposure /
          SEExposure
      )^2
  ]
  
  result[
    ,
    WeakInstrument :=
      !is.na(FStatistic) &
      FStatistic < 10
  ]
  
  result
}


# ------------------------------------------------------------
# Deduplicate metabolite-SNP instruments
# ------------------------------------------------------------

deduplicate_instruments <- function(
    instruments
){
  
  result <- data.table::copy(
    instruments
  )
  
  result[
    ,
    ExposurePForOrdering :=
      data.table::fifelse(
        is.na(PExposure),
        Inf,
        PExposure
      )
  ]
  
  result[
    ,
    FForOrdering :=
      data.table::fifelse(
        is.na(FStatistic),
        -Inf,
        FStatistic
      )
  ]
  
  data.table::setorder(
    result,
    Metabolite,
    SNP,
    ExposurePForOrdering,
    -FForOrdering
  )
  
  result <- unique(
    result,
    by = c(
      "Metabolite",
      "SNP"
    )
  )
  
  result[
    ,
    c(
      "ExposurePForOrdering",
      "FForOrdering"
    ) := NULL
  ]
  
  result
}


# ------------------------------------------------------------
# Finalize InstrumentEvidence table
# ------------------------------------------------------------

build_instrument_evidence <- function(
    harmonized_instruments
){
  
  result <- calculate_missing_instrument_strength(
    harmonized_instruments
  )
  
  result <- deduplicate_instruments(
    result
  )
  
  result[
    ,
    InstrumentStrengthClass :=
      data.table::fcase(
        is.na(FStatistic),
        "Unknown",
        
        FStatistic >= 30,
        "Strong",
        
        FStatistic >= 10,
        "Adequate",
        
        default = "Weak"
      )
  ]
  
  output_columns <- intersect(
    c(
      "Metabolite",
      "StudyAccession",
      "SNP",
      
      "CHR",
      "BP",
      "CHROutcome",
      "BPOutcome",
      
      "EffectAllele",
      "OtherAllele",
      "EffectAlleleFrequency",
      
      "EffectAlleleOutcome",
      "OtherAlleleOutcome",
      "EffectAlleleFrequencyOutcome",
      
      "BetaExposure",
      "SEExposure",
      "PExposure",
      
      "BetaOutcome",
      "SEOutcome",
      "POutcome",
      
      "FStatistic",
      "R2",
      
      "InstrumentStrengthClass",
      "WeakInstrument",
      
      "HarmonizationStatus",
      "Palindromic",
      "Ambiguous",
      "MHC",
      "SteigerDirection",
      
      "SourceFile"
    ),
    names(result)
  )
  
  result <- result[
    ,
    ..output_columns
  ]
  
  data.table::setorder(
    result,
    Metabolite,
    PExposure,
    -FStatistic,
    SNP
  )
  
  result
}


# ------------------------------------------------------------
# Collapse unique character values
# ------------------------------------------------------------

collapse_unique_values <- function(
    x,
    separator = ";"
){
  
  values <- unique(
    as.character(
      x[
        !is.na(x) &
          as.character(x) != ""
      ]
    )
  )
  
  if(length(values) == 0){
    
    return(
      NA_character_
    )
  }
  
  paste(
    values,
    collapse = separator
  )
}


# ------------------------------------------------------------
# Summarize instruments by metabolite
# ------------------------------------------------------------

summarize_instruments_by_metabolite <- function(
    instrument_evidence
){
  
  result <- data.table::copy(
    instrument_evidence
  )
  
  result[
    ,
    LeadOrderingP :=
      data.table::fifelse(
        is.na(PExposure),
        Inf,
        PExposure
      )
  ]
  
  result[
    ,
    LeadOrderingF :=
      data.table::fifelse(
        is.na(FStatistic),
        -Inf,
        FStatistic
      )
  ]
  
  data.table::setorder(
    result,
    Metabolite,
    LeadOrderingP,
    -LeadOrderingF,
    SNP
  )
  
  summary_table <- result[
    ,
    .(
      N_Instruments =
        data.table::uniqueN(SNP),
      
      N_StrongInstruments =
        sum(
          InstrumentStrengthClass == "Strong",
          na.rm = TRUE
        ),
      
      N_AdequateInstruments =
        sum(
          InstrumentStrengthClass == "Adequate",
          na.rm = TRUE
        ),
      
      N_WeakInstruments =
        sum(
          InstrumentStrengthClass == "Weak",
          na.rm = TRUE
        ),
      
      LeadSNP =
        SNP[1],
      
      LeadSNP_CHR =
        CHR[1],
      
      LeadSNP_BP =
        BP[1],
      
      LeadSNP_P =
        PExposure[1],
      
      LeadSNP_F =
        FStatistic[1],
      
      MinimumExposureP =
        suppressWarnings(
          min(
            PExposure,
            na.rm = TRUE
          )
        ),
      
      MedianFStatistic =
        suppressWarnings(
          median(
            FStatistic,
            na.rm = TRUE
          )
        ),
      
      MinimumFStatistic =
        suppressWarnings(
          min(
            FStatistic,
            na.rm = TRUE
          )
        ),
      
      MaximumFStatistic =
        suppressWarnings(
          max(
            FStatistic,
            na.rm = TRUE
          )
        ),
      
      InstrumentSNPs =
        collapse_unique_values(
          SNP
        ),
      
      InstrumentPositions =
        collapse_unique_values(
          paste0(
            CHR,
            ":",
            BP
          )
        ),
      
      AnyMHCInstrument =
        any(
          MHC %in% TRUE,
          na.rm = TRUE
        ),
      
      AnyPalindromicInstrument =
        any(
          Palindromic %in% TRUE,
          na.rm = TRUE
        )
    ),
    by = .(
      Metabolite,
      StudyAccession
    )
  ]
  
  numeric_summary_columns <- c(
    "MinimumExposureP",
    "MedianFStatistic",
    "MinimumFStatistic",
    "MaximumFStatistic"
  )
  
  for(column in numeric_summary_columns){
    
    summary_table[
      !is.finite(
        get(column)
      ),
      (column) := NA_real_
    ]
  }
  
  summary_table[
    ,
    InstrumentEvidenceClass :=
      data.table::fcase(
        N_Instruments == 0,
        "NoInstrument",
        
        !is.na(MinimumFStatistic) &
          MinimumFStatistic >= 10,
        "AllInstrumentsAdequate",
        
        N_WeakInstruments == 0,
        "NoKnownWeakInstrument",
        
        N_WeakInstruments > 0,
        "ContainsWeakInstrument",
        
        default = "UnknownStrength"
      )
  ]
  
  data.table::setorder(
    summary_table,
    Metabolite,
    StudyAccession
  )
  
  summary_table
}


# ------------------------------------------------------------
# Join instrument summary to evidence table
# ------------------------------------------------------------

attach_instrument_summary <- function(
    evidence_table,
    instrument_summary
){
  
  result <- data.table::copy(
    evidence_table
  )
  
  summary_copy <- data.table::copy(
    instrument_summary
  )
  
  # SNP özeti temel olarak metabolit düzeyindedir.
  # Accession farklı kaynaklarda farklı anlam taşıyabileceği için
  # burada yalnızca Metabolite üzerinden bağlanır.
  
  if("StudyAccession" %in% names(summary_copy)){
    
    summary_copy[
      ,
      StudyAccession := NULL
    ]
  }
  
  summary_copy <- unique(
    summary_copy,
    by = "Metabolite"
  )
  
  merge(
    result,
    summary_copy,
    by = "Metabolite",
    all.x = TRUE,
    sort = FALSE
  )
}


# ------------------------------------------------------------
# Build SNP-metabolite-drug evidence table
# ------------------------------------------------------------

# ------------------------------------------------------------
# Build SNP-metabolite-drug-class evidence table
# ------------------------------------------------------------

build_variant_drug_evidence <- function(
    instrument_evidence,
    candidate_network
){
  instrument_copy <- data.table::copy(instrument_evidence)
  candidate_copy <- data.table::copy(candidate_network)

  check_columns(instrument_copy, c("Metabolite", "SNP"))
  check_columns(candidate_copy, c("Metabolite", "DrugCategory"))

  join_columns <- "Metabolite"
  if (
    "StudyAccession" %in% names(instrument_copy) &&
      "study_accession" %in% names(candidate_copy)
  ) {
    candidate_copy[, StudyAccession := as.character(study_accession)]
    join_columns <- c("Metabolite", "StudyAccession")
  }

  variant_drug <- merge(
    instrument_copy,
    candidate_copy,
    by = join_columns,
    allow.cartesian = TRUE,
    sort = FALSE
  )

  preferred_columns <- intersect(
    c(
      "SNP",
      "CHR",
      "BP",
      "EffectAllele",
      "OtherAllele",
      "EffectAlleleFrequency",
      "Metabolite",
      "StudyAccession",
      "BetaExposure",
      "SEExposure",
      "PExposure",
      "BetaOutcome",
      "SEOutcome",
      "POutcome",
      "FStatistic",
      "R2",
      "InstrumentStrengthClass",
      "WeakInstrument",
      "HarmonizationStatus",
      "Palindromic",
      "Ambiguous",
      "MHC",
      "SteigerDirection",
      "NatMedAbbreviation",
      "NatMedExplanation",
      "MetaboliteGroup",
      "reported_trait",
      "MR_beta",
      "MR_se",
      "MR_p",
      "MR_FDR",
      "PP_H4",
      "MR_Direction",
      "DrugCategory",
      "CategoryCode",
      "DrugName",
      "AnatomicalMainGroup",
      "GroupLetter",
      "ATCIncluded",
      "CardiovascularATCGroup",
      "M1_Beta",
      "M1_P",
      "M2_UserN",
      "M2_NonUserN",
      "M2_Beta",
      "M2_SE",
      "M2_P",
      "EffectiveM2Threshold",
      "M2Significant",
      "M1M2Significant",
      "M1M2DirectionConsistent",
      "M1M2Robust",
      "DrugMetaboliteDirection",
      "DirectionalConcordance",
      "DirectionClass",
      "MREvidenceClass",
      "ColocalizationEvidenceClass",
      "PharmacometabolomicEvidenceClass",
      "EvidenceTier",
      "CandidateSelected",
      "CandidateSelectionReason",
      "MetabolitePriorityRank",
      "GlobalPriorityRank"
    ),
    names(variant_drug)
  )

  variant_drug <- variant_drug[, ..preferred_columns]

  order_columns <- intersect(
    c("Metabolite", "MetabolitePriorityRank", "PExposure", "SNP", "DrugCategory"),
    names(variant_drug)
  )
  data.table::setorderv(variant_drug, order_columns)

  variant_drug
}

# ------------------------------------------------------------
# Validate variant evidence outputs
# ------------------------------------------------------------

validate_variant_evidence <- function(
    instrument_evidence,
    instrument_summary,
    variant_drug_evidence
){
  check_columns(instrument_evidence, c("Metabolite", "SNP"))
  check_columns(
    instrument_summary,
    c("Metabolite", "N_Instruments", "LeadSNP", "InstrumentSNPs")
  )
  check_columns(
    variant_drug_evidence,
    c("SNP", "Metabolite", "DrugCategory")
  )

  duplicate_instruments <- instrument_evidence[
    , .N,
    by = .(Metabolite, SNP)
  ][N > 1]

  if (nrow(duplicate_instruments) > 0) {
    stop_framework("InstrumentEvidence contains duplicated metabolite-SNP rows.")
  }

  if (any(instrument_summary$N_Instruments < 1, na.rm = TRUE)) {
    stop_framework("MetaboliteInstrumentSummary contains invalid instrument counts.")
  }

  invisible(TRUE)
}
