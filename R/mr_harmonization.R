# ============================================================
# Framework_v1
# MR harmonization functions -- optimized for cached GWS input
# ============================================================

library(data.table)

read_outcome_gwas <- function(outcome_file){
  check_file_exists(outcome_file)

  outcome <- fread(
    outcome_file,
    select = c(
      "rsid",
      "chromosome",
      "base_pair_location",
      "effect_allele",
      "other_allele",
      "beta",
      "standard_error",
      "effect_allele_frequency",
      "p_value"
    )
  )

  setnames(
    outcome,
    old = c(
      "rsid",
      "chromosome",
      "base_pair_location",
      "effect_allele",
      "other_allele",
      "beta",
      "standard_error",
      "effect_allele_frequency",
      "p_value"
    ),
    new = c(
      "SNP",
      "CHR_outcome",
      "BP_outcome",
      "EA_outcome",
      "OA_outcome",
      "BETA_outcome",
      "SE_outcome",
      "EAF_outcome",
      "P_outcome"
    )
  )

  outcome <- outcome[!is.na(SNP) & SNP != ""]
  outcome[, EA_outcome := toupper(EA_outcome)]
  outcome[, OA_outcome := toupper(OA_outcome)]

  setorder(outcome, P_outcome)
  unique(outcome, by = "SNP")
}

read_exposure_gws_cache <- function(gws_file, snps = NULL){
  check_file_exists(gws_file)

  exposure <- fread(gws_file, showProgress = FALSE)

  required <- c(
    "SNP", "CHR", "BP",
    "effect_allele", "other_allele",
    "beta", "se", "eaf", "P"
  )
  check_columns(exposure, required)

  if (!is.null(snps)) {
    exposure <- exposure[SNP %in% snps]
  }

  exposure <- exposure[!is.na(SNP) & SNP != ""]

  exposure <- exposure[, .(
    SNP = as.character(SNP),
    CHR_exposure = CHR,
    BP_exposure = BP,
    EA_exposure = toupper(as.character(effect_allele)),
    OA_exposure = toupper(as.character(other_allele)),
    BETA_exposure = as.numeric(beta),
    SE_exposure = as.numeric(se),
    EAF_exposure = as.numeric(eaf),
    P_exposure = as.numeric(P)
  )]

  setorder(exposure, P_exposure)
  unique(exposure, by = "SNP")
}

harmonize_alleles <- function(exposure, outcome){
  dat <- merge(
    exposure,
    outcome,
    by = "SNP"
  )

  dat[, harmonisation := fifelse(
    EA_exposure == EA_outcome & OA_exposure == OA_outcome,
    "aligned",
    fifelse(
      EA_exposure == OA_outcome & OA_exposure == EA_outcome,
      "flipped",
      "mismatch"
    )
  )]

  dat[, BETA_outcome_harmonized := BETA_outcome]
  dat[, EAF_outcome_harmonized := EAF_outcome]

  dat[harmonisation == "flipped",
      BETA_outcome_harmonized := -BETA_outcome]

  dat[harmonisation == "flipped" & !is.na(EAF_outcome),
      EAF_outcome_harmonized := 1 - EAF_outcome]

  dat
}

remove_allele_mismatches <- function(dat){
  dat[harmonisation != "mismatch"]
}

remove_ambiguous_palindromic <- function(dat, maf_lower = 0.42, maf_upper = 0.58){
  dat[, palindromic := paste0(EA_exposure, OA_exposure) %in%
        c("AT", "TA", "CG", "GC")]

  dat[, ambiguous_palindromic := palindromic &
        !is.na(EAF_exposure) &
        EAF_exposure > maf_lower &
        EAF_exposure < maf_upper]

  dat[ambiguous_palindromic != TRUE]
}

create_mr_ready_table <- function(dat, metabolite, outcome_name){
  dat[, .(
    SNP,
    beta.exposure = BETA_exposure,
    se.exposure = SE_exposure,
    effect_allele.exposure = EA_exposure,
    other_allele.exposure = OA_exposure,
    eaf.exposure = EAF_exposure,
    pval.exposure = P_exposure,
    beta.outcome = BETA_outcome_harmonized,
    se.outcome = SE_outcome,
    effect_allele.outcome = EA_exposure,
    other_allele.outcome = OA_exposure,
    eaf.outcome = EAF_outcome_harmonized,
    pval.outcome = P_outcome,
    CHR_exposure,
    BP_exposure,
    CHR_outcome,
    BP_outcome,
    harmonisation,
    palindromic,
    exposure = metabolite,
    outcome = outcome_name
  )]
}

prepare_mr_input_for_metabolite <- function(
    metabolite,
    clumped_file,
    gws_instrument_file,
    outcome,
    outcome_name,
    output_dir
){
  create_directory(output_dir)

  check_file_exists(clumped_file)
  check_file_exists(gws_instrument_file)

  clumped <- fread(clumped_file, fill = TRUE, showProgress = FALSE)

  if (!"SNP" %in% names(clumped)) {
    stop_framework(
      paste0(
        "PLINK .clumped output missing SNP column for ",
        metabolite
      )
    )
  }

  snps <- unique(as.character(clumped$SNP))
  snps <- snps[!is.na(snps) & snps != ""]

  exposure <- read_exposure_gws_cache(
    gws_file = gws_instrument_file,
    snps = snps
  )

  dat <- harmonize_alleles(
    exposure = exposure,
    outcome = outcome
  )

  n_rsID_matches <- nrow(dat)
  n_aligned <- nrow(dat[harmonisation == "aligned"])
  n_flipped <- nrow(dat[harmonisation == "flipped"])
  n_mismatch <- nrow(dat[harmonisation == "mismatch"])

  dat <- remove_allele_mismatches(dat)

  n_before_pal <- nrow(dat)
  dat <- remove_ambiguous_palindromic(dat)
  n_pal_removed <- n_before_pal - nrow(dat)

  mr_ready <- create_mr_ready_table(
    dat = dat,
    metabolite = metabolite,
    outcome_name = outcome_name
  )

  out_file <- file.path(
    output_dir,
    paste0(metabolite, "_harmonized_for_MR_rsID.tsv")
  )

  fwrite(mr_ready, out_file, sep = "\t")

  data.table(
    Metabolite = metabolite,
    n_clumped_snps = length(snps),
    n_exposure_snps_found = nrow(exposure),
    n_rsID_matches_with_outcome = n_rsID_matches,
    n_aligned = n_aligned,
    n_flipped = n_flipped,
    n_mismatch = n_mismatch,
    n_palindromic_removed = n_pal_removed,
    n_harmonized_snps = nrow(mr_ready),
    output_file = out_file,
    status = "ok"
  )
}
