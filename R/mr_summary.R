# ============================================================
# Framework_v1
# MR summary functions
# ============================================================

library(data.table)

#--------------------------------------------------------------
# Select IVW results
#--------------------------------------------------------------

extract_ivw_results <- function(mr_results){
  
  if(is.null(mr_results) || nrow(mr_results) == 0){
    
    return(data.table())
    
  }
  
  required_cols <- c(
    "Metabolite",
    "method",
    "nsnp",
    "b",
    "se",
    "pval"
  )
  
  check_columns(mr_results, required_cols)
  
  ivw <- mr_results[
    method == "Inverse variance weighted"
  ]
  
  if(nrow(ivw) == 0){
    
    return(data.table())
    
  }
  
  ivw <- copy(ivw)
  
  ivw[, MR_beta := b]
  ivw[, MR_se := se]
  ivw[, MR_p := pval]
  
  if(!"OR" %in% names(ivw)){
    ivw[, OR := exp(MR_beta)]
  }
  
  if(!"OR_lci95" %in% names(ivw)){
    ivw[, OR_lci95 := exp(MR_beta - 1.96 * MR_se)]
  }
  
  if(!"OR_uci95" %in% names(ivw)){
    ivw[, OR_uci95 := exp(MR_beta + 1.96 * MR_se)]
  }
  
  ivw[, .(
    Metabolite,
    nsnp,
    MR_beta,
    MR_se,
    MR_p,
    OR,
    OR_lci95,
    OR_uci95
  )]
  
}

#--------------------------------------------------------------
# Add multiple-testing correction
#--------------------------------------------------------------

add_mr_fdr <- function(ivw_results, method = "BH"){
  
  result <- copy(ivw_results)
  
  if(nrow(result) == 0){
    return(result)
  }
  
  result[, MR_FDR := p.adjust(
    MR_p,
    method = method
  )]
  
  result
  
}

#--------------------------------------------------------------
# Extract IVW heterogeneity results
#--------------------------------------------------------------

extract_ivw_heterogeneity <- function(heterogeneity_results){
  
  if(
    is.null(heterogeneity_results) ||
    nrow(heterogeneity_results) == 0
  ){
    
    return(data.table())
    
  }
  
  required_cols <- c(
    "Metabolite",
    "method",
    "Q",
    "Q_df",
    "Q_pval"
  )
  
  if(!all(required_cols %in% names(heterogeneity_results))){
    return(data.table())
  }
  
  result <- heterogeneity_results[
    method == "Inverse variance weighted"
  ]
  
  result[, .(
    Metabolite,
    Heterogeneity_Q = Q,
    Heterogeneity_df = Q_df,
    Heterogeneity_p = Q_pval
  )]
  
}

#--------------------------------------------------------------
# Extract Egger intercept results
#--------------------------------------------------------------

extract_egger_intercept <- function(pleiotropy_results){
  
  if(
    is.null(pleiotropy_results) ||
    nrow(pleiotropy_results) == 0
  ){
    
    return(data.table())
    
  }
  
  required_cols <- c(
    "Metabolite",
    "egger_intercept",
    "se",
    "pval"
  )
  
  if(!all(required_cols %in% names(pleiotropy_results))){
    return(data.table())
  }
  
  pleiotropy_results[, .(
    Metabolite,
    EggerIntercept = egger_intercept,
    EggerIntercept_SE = se,
    EggerIntercept_p = pval
  )]
  
}

#--------------------------------------------------------------
# Assign direction of MR effect
#--------------------------------------------------------------

assign_mr_direction <- function(mr_beta){
  
  direction <- rep(
    NA_character_,
    length(mr_beta)
  )
  
  direction[
    !is.na(mr_beta) &
      mr_beta > 0
  ] <- "Risk_increasing"
  
  direction[
    !is.na(mr_beta) &
      mr_beta < 0
  ] <- "Risk_decreasing"
  
  direction[
    !is.na(mr_beta) &
      mr_beta == 0
  ] <- "Neutral"
  
  direction
  
}

#--------------------------------------------------------------
# Assign MR priority category
#--------------------------------------------------------------

assign_mr_priority <- function(
    mr_fdr,
    mr_p,
    n_strong,
    egger_p,
    fdr_threshold = 0.05,
    nominal_threshold = 0.05,
    minimum_strong_instruments = 3
){
  
  result <- rep(
    "low_priority",
    length(mr_p)
  )
  
  high <- (
    !is.na(mr_fdr) &
      mr_fdr < fdr_threshold &
      !is.na(n_strong) &
      n_strong >= minimum_strong_instruments &
      (
        is.na(egger_p) |
          egger_p > 0.05
      )
  )
  
  nominal <- (
    !high &
      !is.na(mr_p) &
      mr_p < nominal_threshold &
      !is.na(n_strong) &
      n_strong >= minimum_strong_instruments
  )
  
  result[high] <- "high_priority"
  result[nominal] <- "nominal_priority"
  
  result
  
}

#--------------------------------------------------------------
# Build metabolite-level MR summary
#--------------------------------------------------------------

build_mr_summary <- function(
    mr_results,
    heterogeneity_results = NULL,
    pleiotropy_results = NULL,
    instrument_summary = NULL,
    mhc_summary = NULL,
    presso_summary = NULL,
    steiger_summary = NULL
){
  
  ivw <- extract_ivw_results(
    mr_results
  )
  
  if(nrow(ivw) == 0){
    
    stop_framework(
      "No IVW results available for MR summary."
    )
    
  }
  
  summary <- add_mr_fdr(ivw)
  
  heterogeneity <- extract_ivw_heterogeneity(
    heterogeneity_results
  )
  
  if(nrow(heterogeneity) > 0){
    
    summary <- merge(
      summary,
      heterogeneity,
      by = "Metabolite",
      all.x = TRUE
    )
    
  }
  
  pleiotropy <- extract_egger_intercept(
    pleiotropy_results
  )
  
  if(nrow(pleiotropy) > 0){
    
    summary <- merge(
      summary,
      pleiotropy,
      by = "Metabolite",
      all.x = TRUE
    )
    
  }
  
  if(
    !is.null(instrument_summary) &&
    nrow(instrument_summary) > 0
  ){
    
    summary <- merge(
      summary,
      instrument_summary,
      by = "Metabolite",
      all.x = TRUE
    )
    
  }
  
  if(
    !is.null(mhc_summary) &&
    nrow(mhc_summary) > 0
  ){
    
    summary <- merge(
      summary,
      mhc_summary,
      by = "Metabolite",
      all.x = TRUE,
      suffixes = c("", "_MHC")
    )
    
  }
  
  if(
    !is.null(presso_summary) &&
    nrow(presso_summary) > 0
  ){
    
    presso_cols <- intersect(
      c(
        "Metabolite",
        "MRPRESSO_GlobalP",
        "MRPRESSO_NOutliers"
      ),
      names(presso_summary)
    )
    
    if(length(presso_cols) > 1){
      
      summary <- merge(
        summary,
        presso_summary[, ..presso_cols],
        by = "Metabolite",
        all.x = TRUE
      )
      
    }
    
  }
  
  if(
    !is.null(steiger_summary) &&
    nrow(steiger_summary) > 0
  ){
    
    steiger_cols <- intersect(
      c(
        "Metabolite",
        "SteigerStatus",
        "SteigerReason"
      ),
      names(steiger_summary)
    )
    
    if(length(steiger_cols) > 1){
      
      summary <- merge(
        summary,
        steiger_summary[, ..steiger_cols],
        by = "Metabolite",
        all.x = TRUE
      )
      
    }
    
  }
  
  summary[, Direction := assign_mr_direction(
    MR_beta
  )]
  
  if(
    all(
      c(
        "n_F_gt_threshold",
        "EggerIntercept_p"
      ) %in% names(summary)
    )
  ){
    
    summary[, MR_priority := assign_mr_priority(
      mr_fdr = MR_FDR,
      mr_p = MR_p,
      n_strong = n_F_gt_threshold,
      egger_p = EggerIntercept_p
    )]
    
  } else {
    
    summary[, MR_priority := NA_character_]
    
  }
  
  setorder(
    summary,
    MR_FDR,
    MR_p
  )
  
  summary
  
}

#--------------------------------------------------------------
# Validate MR summary structure
#--------------------------------------------------------------

validate_mr_summary <- function(mr_summary){
  
  required_cols <- c(
    "Metabolite",
    "nsnp",
    "MR_beta",
    "MR_se",
    "MR_p",
    "MR_FDR",
    "Direction"
  )
  
  check_columns(
    mr_summary,
    required_cols
  )
  
  if(anyDuplicated(mr_summary$Metabolite)){
    
    stop_framework(
      "MR summary contains duplicated metabolites."
    )
    
  }
  
  invisible(TRUE)
  
}
