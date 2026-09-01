# ============================================================
# Framework_v1
# MR analysis functions
# ============================================================

library(data.table)
library(TwoSampleMR)

#--------------------------------------------------------------
# Convert MR-ready table to TwoSampleMR format
#--------------------------------------------------------------

prepare_twosamplemr_input <- function(
    mr_ready,
    metabolite,
    outcome_name
){
  
  required_cols <- c(
    "SNP",
    "beta.exposure",
    "se.exposure",
    "effect_allele.exposure",
    "other_allele.exposure",
    "eaf.exposure",
    "pval.exposure",
    "beta.outcome",
    "se.outcome",
    "effect_allele.outcome",
    "other_allele.outcome",
    "eaf.outcome",
    "pval.outcome"
  )
  
  check_columns(mr_ready, required_cols)
  
  base_cols <- c(
    "SNP",
    "beta.exposure",
    "se.exposure",
    "effect_allele.exposure",
    "other_allele.exposure",
    "eaf.exposure",
    "pval.exposure",
    "beta.outcome",
    "se.outcome",
    "effect_allele.outcome",
    "other_allele.outcome",
    "eaf.outcome",
    "pval.outcome"
  )
  
  optional_cols <- intersect(
    c(
      "CHR_exposure",
      "BP_exposure",
      "CHR_outcome",
      "BP_outcome",
      "harmonisation",
      "palindromic"
    ),
    names(mr_ready)
  )
  
  harm <- mr_ready[, c(
    base_cols,
    optional_cols
  ), with = FALSE]
  
  harm[, exposure := metabolite]
  harm[, outcome := outcome_name]
  
  harm[, id.exposure := metabolite]
  harm[, id.outcome := outcome_name]
  
  harm[, mr_keep := TRUE]
  
  harm[, pval.exposure := as.numeric(pval.exposure)]
  harm[, pval.outcome := as.numeric(pval.outcome)]
  
  harm
  
}

#--------------------------------------------------------------
# Run MR methods
#--------------------------------------------------------------

run_mr_methods <- function(
    harm,
    method_list = c(
      "mr_ivw",
      "mr_egger_regression",
      "mr_weighted_median",
      "mr_weighted_mode",
      "mr_simple_mode"
    )
){
  
  if(nrow(harm) < 3){
    
    return(data.table(
      status = "skipped",
      reason = "fewer_than_3_snps"
    ))
    
  }
  
  result <- tryCatch(
    {
      as.data.table(
        TwoSampleMR::mr(
          harm,
          method_list = method_list
        )
      )
    },
    error = function(e){
      
      data.table(
        status = "error",
        error = e$message
      )
      
    }
  )
  
  result
  
}

#--------------------------------------------------------------
# Add odds ratios for binary outcome
#--------------------------------------------------------------

add_odds_ratios <- function(mr_results){
  
  required_cols <- c("b", "se")
  
  if(!all(required_cols %in% names(mr_results))){
    return(mr_results)
  }
  
  mr_results[, OR := exp(b)]
  mr_results[, OR_lci95 := exp(b - 1.96 * se)]
  mr_results[, OR_uci95 := exp(b + 1.96 * se)]
  
  mr_results
  
}

#--------------------------------------------------------------
# Run heterogeneity analysis
#--------------------------------------------------------------

run_mr_heterogeneity <- function(harm){
  
  if(nrow(harm) < 3){
    
    return(data.table(
      status = "skipped",
      reason = "fewer_than_3_snps"
    ))
    
  }
  
  tryCatch(
    {
      as.data.table(
        TwoSampleMR::mr_heterogeneity(harm)
      )
    },
    error = function(e){
      
      data.table(
        status = "error",
        error = e$message
      )
      
    }
  )
  
}

#--------------------------------------------------------------
# Run Egger intercept pleiotropy test
#--------------------------------------------------------------

run_mr_pleiotropy <- function(harm){
  
  if(nrow(harm) < 3){
    
    return(data.table(
      status = "skipped",
      reason = "fewer_than_3_snps"
    ))
    
  }
  
  tryCatch(
    {
      as.data.table(
        TwoSampleMR::mr_pleiotropy_test(harm)
      )
    },
    error = function(e){
      
      data.table(
        status = "error",
        error = e$message
      )
      
    }
  )
  
}

#--------------------------------------------------------------
# Run complete MR analysis for one metabolite
#--------------------------------------------------------------

run_mr_for_metabolite <- function(
    mr_ready,
    metabolite,
    outcome_name,
    method_list = c(
      "mr_ivw",
      "mr_egger_regression",
      "mr_weighted_median",
      "mr_weighted_mode",
      "mr_simple_mode"
    )
){
  
  harm <- prepare_twosamplemr_input(
    mr_ready = mr_ready,
    metabolite = metabolite,
    outcome_name = outcome_name
  )
  
  mr_results <- run_mr_methods(
    harm = harm,
    method_list = method_list
  )
  
  mr_results <- add_odds_ratios(mr_results)
  
  if("status" %in% names(mr_results)){
    mr_results[, Metabolite := metabolite]
  } else {
    mr_results[, Metabolite := metabolite]
  }
  
  heterogeneity <- run_mr_heterogeneity(harm)
  heterogeneity[, Metabolite := metabolite]
  
  pleiotropy <- run_mr_pleiotropy(harm)
  pleiotropy[, Metabolite := metabolite]
  
  list(
    input = harm,
    mr = mr_results,
    heterogeneity = heterogeneity,
    pleiotropy = pleiotropy
  )
  
}
