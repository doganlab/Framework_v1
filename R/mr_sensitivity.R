# ============================================================
# Framework_v1
# MR sensitivity and quality-control functions
# ============================================================

library(data.table)

#--------------------------------------------------------------
# Calculate SNP-level F-statistics
#--------------------------------------------------------------

calculate_f_statistics <- function(harm){
  
  required_cols <- c(
    "beta.exposure",
    "se.exposure"
  )
  
  check_columns(harm, required_cols)
  
  result <- copy(harm)
  
  result[, F_stat := (beta.exposure / se.exposure)^2]
  
  result[, InstrumentStrength := fifelse(
    is.na(F_stat),
    NA_character_,
    fifelse(
      F_stat > 10,
      "strong",
      "weak"
    )
  )]
  
  result
  
}

#--------------------------------------------------------------
# Summarize instrument strength
#--------------------------------------------------------------

summarize_instrument_strength <- function(
    harm,
    metabolite,
    f_threshold = 10
){
  
  dat <- calculate_f_statistics(harm)
  
  valid_f <- dat$F_stat[
    is.finite(dat$F_stat)
  ]
  
  data.table(
    Metabolite = metabolite,
    n_snps_total = nrow(dat),
    n_F_gt_threshold = sum(
      dat$F_stat > f_threshold,
      na.rm = TRUE
    ),
    n_F_le_threshold = sum(
      dat$F_stat <= f_threshold,
      na.rm = TRUE
    ),
    mean_F = if(length(valid_f) > 0) {
      mean(valid_f)
    } else {
      NA_real_
    },
    median_F = if(length(valid_f) > 0) {
      median(valid_f)
    } else {
      NA_real_
    },
    min_F = if(length(valid_f) > 0) {
      min(valid_f)
    } else {
      NA_real_
    },
    max_F = if(length(valid_f) > 0) {
      max(valid_f)
    } else {
      NA_real_
    }
  )
  
}

#--------------------------------------------------------------
# Retain strong instruments
#--------------------------------------------------------------

filter_strong_instruments <- function(
    harm,
    f_threshold = 10
){
  
  dat <- calculate_f_statistics(harm)
  
  dat[
    !is.na(F_stat) &
      F_stat > f_threshold
  ]
  
}

#--------------------------------------------------------------
# Check instruments located in the MHC region
#--------------------------------------------------------------

summarize_mhc_instruments <- function(
    harm,
    metabolite,
    mhc_chr = 6,
    mhc_start = 25000000,
    mhc_end = 34000000
){
  
  required_cols <- c(
    "CHR_exposure",
    "BP_exposure"
  )
  
  check_columns(harm, required_cols)
  
  chr_values <- as.character(harm$CHR_exposure)
  chr_values <- gsub(
    "^chr",
    "",
    chr_values,
    ignore.case = TRUE
  )
  
  bp_values <- as.numeric(harm$BP_exposure)
  
  is_mhc <- (
    chr_values == as.character(mhc_chr) &
      !is.na(bp_values) &
      bp_values >= mhc_start &
      bp_values <= mhc_end
  )
  
  n_snps <- nrow(harm)
  n_mhc <- sum(is_mhc, na.rm = TRUE)
  
  data.table(
    Metabolite = metabolite,
    n_snps = n_snps,
    n_chr6 = sum(
      chr_values == as.character(mhc_chr),
      na.rm = TRUE
    ),
    n_mhc = n_mhc,
    mhc_percent = if(n_snps > 0) {
      100 * n_mhc / n_snps
    } else {
      NA_real_
    }
  )
  
}

#--------------------------------------------------------------
# Run leave-one-out analysis
#--------------------------------------------------------------

run_leave_one_out <- function(
    harm,
    metabolite,
    minimum_snps = 3
){
  
  if(nrow(harm) < minimum_snps){
    
    return(data.table(
      Metabolite = metabolite,
      status = "skipped",
      reason = paste0(
        "fewer_than_",
        minimum_snps,
        "_snps"
      )
    ))
    
  }
  
  result <- tryCatch(
    {
      x <- TwoSampleMR::mr_leaveoneout(harm)
      x <- as.data.table(x)
      
      x[, Metabolite := metabolite]
      x[, status := "ok"]
      
      x
    },
    error = function(e){
      
      data.table(
        Metabolite = metabolite,
        status = "error",
        error = e$message
      )
      
    }
  )
  
  result
  
}

#--------------------------------------------------------------
# Run MR-PRESSO
#--------------------------------------------------------------

run_mr_presso <- function(
    harm,
    metabolite,
    minimum_snps = 4,
    nb_distribution = 1000,
    significance_threshold = 0.05
){
  
  if(nrow(harm) < minimum_snps){
    
    return(data.table(
      Metabolite = metabolite,
      status = "skipped",
      reason = paste0(
        "fewer_than_",
        minimum_snps,
        "_snps"
      ),
      MRPRESSO_GlobalP = NA_real_,
      MRPRESSO_NOutliers = NA_integer_
    ))
    
  }
  
  if(!requireNamespace("MRPRESSO", quietly = TRUE)){
    
    return(data.table(
      Metabolite = metabolite,
      status = "skipped",
      reason = "MRPRESSO_package_not_installed",
      MRPRESSO_GlobalP = NA_real_,
      MRPRESSO_NOutliers = NA_integer_
    ))
    
  }
  
  result <- tryCatch(
    {
      
      presso <- MRPRESSO::mr_presso(
        BetaOutcome = "beta.outcome",
        BetaExposure = "beta.exposure",
        SdOutcome = "se.outcome",
        SdExposure = "se.exposure",
        OUTLIERtest = TRUE,
        DISTORTIONtest = TRUE,
        data = as.data.frame(harm),
        NbDistribution = nb_distribution,
        SignifThreshold = significance_threshold
      )
      
      global_test <- presso[["MR-PRESSO results"]][["Global Test"]]
      
      global_p <- NA_real_
      
      if (!is.null(global_test)) {
        
        if (is.data.frame(global_test) || is.list(global_test)) {
          
          p_name <- grep(
            "Pvalue|P-value|pvalue|P.value|P-value",
            names(global_test),
            value = TRUE,
            ignore.case = TRUE
          )[1]
          
          if (!is.na(p_name) && length(p_name) > 0) {
            global_p <- suppressWarnings(
              as.numeric(global_test[[p_name]][1])
            )
          }
          
        } else if (is.numeric(global_test)) {
          
          global_p <- as.numeric(global_test[1])
        }
      }
      
      outlier_test <- presso[["MR-PRESSO results"]][["Outlier Test"]]
      
      n_outliers <- 0L
      
      if(
        !is.null(outlier_test) &&
        is.data.frame(outlier_test) &&
        nrow(outlier_test) > 0
      ){
        
        p_name <- grep(
          "Pvalue|P-value|pvalue",
          names(outlier_test),
          value = TRUE
        )[1]
        
        if(!is.na(p_name)){
          
          outlier_p <- as.numeric(
            outlier_test[[p_name]]
          )
          
          n_outliers <- sum(
            outlier_p <
              significance_threshold /
              nrow(outlier_test),
            na.rm = TRUE
          )
          
        }
        
      }
      
      data.table(
        Metabolite = metabolite,
        status = "ok",
        MRPRESSO_GlobalP = global_p,
        MRPRESSO_NOutliers = n_outliers
      )
      
    },
    error = function(e){
      
      data.table(
        Metabolite = metabolite,
        status = "error",
        error = e$message,
        MRPRESSO_GlobalP = NA_real_,
        MRPRESSO_NOutliers = NA_integer_
      )
      
    }
  )
  
  result
  
}

#--------------------------------------------------------------
# Record Steiger directionality status
#--------------------------------------------------------------

assess_steiger_availability <- function(
    harm,
    metabolite
){
  
  required_cols <- c(
    "samplesize.exposure",
    "samplesize.outcome"
  )
  
  if(!all(required_cols %in% names(harm))){
    
    return(data.table(
      Metabolite = metabolite,
      SteigerStatus = "not_run",
      SteigerReason =
        "sample_size_columns_not_available"
    ))
    
  }
  
  if(
    all(is.na(harm$samplesize.exposure)) ||
    all(is.na(harm$samplesize.outcome))
  ){
    
    return(data.table(
      Metabolite = metabolite,
      SteigerStatus = "not_run",
      SteigerReason =
        "sample_sizes_are_missing"
    ))
    
  }
  
  data.table(
    Metabolite = metabolite,
    SteigerStatus = "available",
    SteigerReason = NA_character_
  )
  
}

#--------------------------------------------------------------
# Run sensitivity analyses for one metabolite
#--------------------------------------------------------------

run_sensitivity_for_metabolite <- function(
    harm,
    metabolite,
    f_threshold = 10,
    run_presso = TRUE,
    presso_nb_distribution = 1000
){
  
  f_stat_data <- calculate_f_statistics(harm)
  
  instrument_summary <- summarize_instrument_strength(
    harm = harm,
    metabolite = metabolite,
    f_threshold = f_threshold
  )
  
  mhc_summary <- summarize_mhc_instruments(
    harm = harm,
    metabolite = metabolite
  )
  
  leave_one_out <- run_leave_one_out(
    harm = harm,
    metabolite = metabolite
  )
  
  if(run_presso){
    
    presso <- run_mr_presso(
      harm = harm,
      metabolite = metabolite,
      nb_distribution =
        presso_nb_distribution
    )
    
  } else {
    
    presso <- data.table(
      Metabolite = metabolite,
      status = "skipped",
      reason = "disabled_by_configuration",
      MRPRESSO_GlobalP = NA_real_,
      MRPRESSO_NOutliers = NA_integer_
    )
    
  }
  
  steiger <- assess_steiger_availability(
    harm = harm,
    metabolite = metabolite
  )
  
  list(
    instrument_data = f_stat_data,
    instrument_summary = instrument_summary,
    mhc_summary = mhc_summary,
    leave_one_out = leave_one_out,
    mr_presso = presso,
    steiger = steiger
  )
  
}

