#!/usr/bin/env Rscript

# ============================================================
# Framework_v1
# 01_MR_pipeline.R -- NatMed / cached-GWS optimized version
# ============================================================

source(file.path(
  dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), mustWork = TRUE)),
  "..", "R", "bootstrap.R"
))

module_files <- c(
  file.path(PROJECT_ROOT, "R", "io.R"),
  file.path(PROJECT_ROOT, "R", "utils.R"),
  file.path(PROJECT_ROOT, "R", "mr_input.R"),
  file.path(PROJECT_ROOT, "R", "mr_clumping.R"),
  file.path(PROJECT_ROOT, "R", "mr_harmonization.R"),
  file.path(PROJECT_ROOT, "R", "mr_analysis.R")
)

for(module_file in module_files){
  if(!file.exists(module_file)){
    stop(paste0("Required module not found: ", module_file), call. = FALSE)
  }
  source(module_file)
}

required_packages <- c("data.table", "TwoSampleMR")
if (isTRUE(MR_RUN_PRESSO)) {
  required_packages <- c(required_packages, "MRPRESSO")
}

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if(length(missing_packages) > 0){
  stop(
    paste0(
      "Missing R packages:\n",
      paste(missing_packages, collapse = "\n")
    ),
    call. = FALSE
  )
}

check_file_exists(OUTCOME_GWAS_FILE)
check_file_exists(CANDIDATE_METABOLITE_FILE)
check_dir_exists(MGWAS_DIR)

ld_required_files <- paste0(
  LD_REFERENCE_PREFIX,
  c(".bed", ".bim", ".fam")
)
missing_ld_files <- ld_required_files[!file.exists(ld_required_files)]

if(length(missing_ld_files) > 0){
  stop_framework(
    paste0(
      "Missing LD reference files:\n",
      paste(missing_ld_files, collapse = "\n")
    )
  )
}

output_directories <- c(
  MR_GWS_DIR,
  MR_CLUMP_DIR,
  MR_HARMONIZED_DIR,
  MR_RESULT_DIR,
  MR_SENSITIVITY_DIR,
  DIR_PROCESSED,
  DIR_LOGS
)

invisible(lapply(output_directories, create_directory))

pipeline_start <- Sys.time()

log_file <- file.path(
  DIR_LOGS,
  paste0(
    "01_MR_pipeline_",
    format(pipeline_start, "%Y%m%d_%H%M%S"),
    ".log"
  )
)

log_pipeline <- function(message){
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- paste0("[", timestamp, "] ", message)
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
}

log_pipeline("Framework_v1 MR pipeline started.")
log_pipeline(paste0("Project root: ", PROJECT_ROOT))
log_pipeline(paste0("Outcome: ", OUTCOME_NAME))

# ------------------------------------------------------------
# 1. Candidate metabolites
# ------------------------------------------------------------

candidates <- load_candidate_metabolites(
  CANDIDATE_METABOLITE_FILE
)
candidates <- data.table::as.data.table(candidates)

log_pipeline(
  paste0("Candidate metabolites loaded: ", nrow(candidates))
)

candidates <- locate_candidate_mgwas_files(
  candidates = candidates,
  mgwas_dir = MGWAS_DIR
)

n_missing_mgwas <- sum(
  is.na(candidates$mGWAS_file) |
    candidates$mGWAS_file == ""
)

log_pipeline(
  paste0("mGWAS files missing: ", n_missing_mgwas)
)

write_tsv_safe(
  candidates,
  file.path(MR_GWS_DIR, "candidate_mGWAS_manifest.tsv")
)

# ------------------------------------------------------------
# 2. GWS extraction with resume-safe cache
# ------------------------------------------------------------

log_pipeline("Preparing genome-wide significant instruments.")

gws_summary <- prepare_all_clumping_inputs(
  candidates = candidates,
  output_dir = MR_GWS_DIR,
  p_threshold = MR_P_THRESHOLD
)

write_tsv_safe(
  gws_summary,
  file.path(MR_GWS_DIR, "GWS_instrument_preparation_summary.tsv")
)

n_clump_ready <- sum(
  gws_summary$status == "ok",
  na.rm = TRUE
)

log_pipeline(
  paste0("Metabolites ready for clumping: ", n_clump_ready)
)

if(n_clump_ready == 0){
  stop_framework("No metabolites produced valid clumping input files.")
}

# ------------------------------------------------------------
# 3. PLINK clumping
# ------------------------------------------------------------

log_pipeline("Running PLINK LD clumping.")

clumping_manifest <- run_plink_clump_batch(
  input_dir = MR_GWS_DIR,
  output_dir = MR_CLUMP_DIR,
  plink_bin = PLINK_BIN,
  ld_ref_prefix = LD_REFERENCE_PREFIX,
  pattern = "_GWS_for_clumping.tsv$",
  clump_p1 = MR_CLUMP_P1,
  clump_p2 = MR_CLUMP_P2,
  clump_r2 = MR_CLUMP_R2,
  clump_kb = MR_CLUMP_KB
)

# Attach GWS cache path so harmonization never re-reads raw mGWAS.
clumping_manifest <- merge(
  clumping_manifest,
  gws_summary[, .(
    Metabolite,
    gws_instrument_file
  )],
  by = "Metabolite",
  all.x = TRUE
)

write_tsv_safe(
  clumping_manifest,
  file.path(MR_CLUMP_DIR, "clumping_manifest.tsv")
)

clumped_summary <- summarize_clumped_instruments(MR_CLUMP_DIR)

write_tsv_safe(
  clumped_summary,
  file.path(MR_CLUMP_DIR, "clumped_instrument_counts.tsv")
)

log_pipeline(
  paste0("Clumping completed for ", nrow(clumped_summary), " metabolites.")
)

# ------------------------------------------------------------
# 4. CAD outcome
# ------------------------------------------------------------

log_pipeline("Reading outcome GWAS.")

outcome <- read_outcome_gwas(
  OUTCOME_GWAS_FILE
)

log_pipeline(
  paste0(
    "Outcome variants retained after standardization: ",
    nrow(outcome)
  )
)

# ------------------------------------------------------------
# 5. Harmonization -- resume safe
# ------------------------------------------------------------

log_pipeline("Preparing harmonized MR inputs.")

harmonization_summaries <- list()

for(i in seq_len(nrow(clumping_manifest))){

  metabolite <- clumping_manifest$Metabolite[i]
  clumped_file <- clumping_manifest$ClumpedFile[i]
  gws_file <- clumping_manifest$gws_instrument_file[i]

  expected_output <- file.path(
    MR_HARMONIZED_DIR,
    paste0(metabolite, "_harmonized_for_MR_rsID.tsv")
  )

  # Resume-safe: if a valid harmonized file already exists, reuse it.
  if(file.exists(expected_output) && file.info(expected_output)$size > 20){

    existing <- tryCatch(
      fread(expected_output, nrows = MR_MINIMUM_SNPS, showProgress = FALSE),
      error = function(e) NULL
    )

    if(!is.null(existing) && nrow(existing) >= MR_MINIMUM_SNPS){

      log_pipeline(
        paste0("Reusing cached harmonized file: ", metabolite)
      )

      harmonization_summaries[[metabolite]] <- data.table(
        Metabolite = metabolite,
        n_clumped_snps = NA_integer_,
        n_exposure_snps_found = NA_integer_,
        n_rsID_matches_with_outcome = NA_integer_,
        n_aligned = NA_integer_,
        n_flipped = NA_integer_,
        n_mismatch = NA_integer_,
        n_palindromic_removed = NA_integer_,
        n_harmonized_snps = nrow(fread(expected_output, showProgress = FALSE)),
        output_file = expected_output,
        status = "reused"
      )

      next
    }
  }

  log_pipeline(paste0("Harmonizing: ", metabolite))

  harmonization_summary <- tryCatch(
    {
      prepare_mr_input_for_metabolite(
        metabolite = metabolite,
        clumped_file = clumped_file,
        gws_instrument_file = gws_file,
        outcome = outcome,
        outcome_name = OUTCOME_NAME,
        output_dir = MR_HARMONIZED_DIR
      )
    },
    error = function(e){
      warn_framework(
        paste0(
          "Harmonization failed for ",
          metabolite,
          ": ",
          e$message
        )
      )

      data.table(
        Metabolite = metabolite,
        status = "error",
        error = e$message
      )
    }
  )

  harmonization_summaries[[metabolite]] <- harmonization_summary
}

harmonization_summary_dt <- rbindlist(
  harmonization_summaries,
  fill = TRUE,
  use.names = TRUE
)

write_tsv_safe(
  harmonization_summary_dt,
  file.path(MR_HARMONIZED_DIR, "MR_harmonization_summary.tsv")
)

successful_harmonization <- harmonization_summary_dt[
  !is.na(output_file) &
    file.exists(output_file) &
    n_harmonized_snps >= MR_MINIMUM_SNPS
]

log_pipeline(
  paste0(
    "Metabolites with MR-ready harmonized data: ",
    nrow(successful_harmonization)
  )
)

if(nrow(successful_harmonization) == 0){
  stop_framework("No metabolites had sufficient harmonized instruments for MR.")
}

# Full outcome no longer needed.
rm(outcome)
invisible(gc())

# ------------------------------------------------------------
# 6. MR -- resume safe
# ------------------------------------------------------------

all_mr <- list()
all_het <- list()
all_pleio <- list()

for(i in seq_len(nrow(successful_harmonization))){

  metabolite <- successful_harmonization$Metabolite[i]
  input_file <- successful_harmonization$output_file[i]

  mr_file <- file.path(
    MR_RESULT_DIR,
    paste0(metabolite, "_MR_results.tsv")
  )
  het_file <- file.path(
    MR_SENSITIVITY_DIR,
    paste0(metabolite, "_heterogeneity.tsv")
  )
  pleio_file <- file.path(
    MR_SENSITIVITY_DIR,
    paste0(metabolite, "_Egger_intercept.tsv")
  )

  # Resume-safe: reuse finished MR result.
  if(file.exists(mr_file) && file.info(mr_file)$size > 20){

    log_pipeline(
      paste0("Reusing cached MR result: ", metabolite)
    )

    all_mr[[metabolite]] <- fread(mr_file, showProgress = FALSE)

    if(file.exists(het_file)){
      all_het[[metabolite]] <- fread(het_file, showProgress = FALSE)
    }

    if(file.exists(pleio_file)){
      all_pleio[[metabolite]] <- fread(pleio_file, showProgress = FALSE)
    }

    next
  }

  log_pipeline(paste0("Running MR: ", metabolite))

  mr_ready <- fread(
    input_file,
    showProgress = FALSE
  )

  analysis <- run_mr_for_metabolite(
    mr_ready = mr_ready,
    metabolite = metabolite,
    outcome_name = OUTCOME_NAME
  )

  mr_res <- as.data.table(analysis$mr)
  het_res <- as.data.table(analysis$heterogeneity)
  pleio_res <- as.data.table(analysis$pleiotropy)

  fwrite(mr_res, mr_file, sep = "\t")
  fwrite(het_res, het_file, sep = "\t")
  fwrite(pleio_res, pleio_file, sep = "\t")

  all_mr[[metabolite]] <- mr_res
  all_het[[metabolite]] <- het_res
  all_pleio[[metabolite]] <- pleio_res
}

mr_dt <- rbindlist(all_mr, fill = TRUE)
het_dt <- rbindlist(all_het, fill = TRUE)
pleio_dt <- rbindlist(all_pleio, fill = TRUE)

write_tsv_safe(
  mr_dt,
  file.path(MR_RESULT_DIR, "TwoSampleMR_all_metabolites_results.tsv")
)

write_tsv_safe(
  het_dt,
  file.path(MR_SENSITIVITY_DIR, "TwoSampleMR_heterogeneity_results.tsv")
)

write_tsv_safe(
  pleio_dt,
  file.path(MR_SENSITIVITY_DIR, "TwoSampleMR_pleiotropy_Egger_intercept.tsv")
)

# ------------------------------------------------------------
# 7. Primary MR summary
# ------------------------------------------------------------

ivw <- mr_dt[
  method == "Inverse variance weighted"
]

if(nrow(ivw) == 0L){
  stop_framework("No IVW MR results were produced.")
}

ivw[, Metabolite := exposure]

ivw[, MR_FDR := p.adjust(
  pval,
  method = "BH"
)]

ivw[, `:=`(
  MR_beta = b,
  MR_se = se,
  MR_p = pval
)]

if(!"OR" %in% names(ivw)){
  ivw[, `:=`(
    OR = exp(MR_beta),
    OR_lci95 = exp(MR_beta - 1.96 * MR_se),
    OR_uci95 = exp(MR_beta + 1.96 * MR_se)
  )]
}

# Add PP.H4 from candidate manifest.
candidate_meta <- unique(
  candidates[, .(
    Metabolite,
    PP_H4
  )]
)

ivw <- merge(
  ivw,
  candidate_meta,
  by = "Metabolite",
  all.x = TRUE
)

setorder(ivw, MR_p)

write_tsv_safe(
  ivw,
  MR_SUMMARY_FILE
)

log_pipeline(
  paste0(
    "MR pipeline completed. IVW metabolites: ",
    nrow(ivw)
  )
)

log_pipeline(
  paste0(
    "MR summary written to: ",
    MR_SUMMARY_FILE
  )
)

cat("\n============================================================\n")
cat("Framework_v1 MR pipeline completed successfully.\n")
cat("============================================================\n")
cat("Candidates:          ", nrow(candidates), "\n")
cat("Clump-ready:         ", n_clump_ready, "\n")
cat("MR-ready:            ", nrow(successful_harmonization), "\n")
cat("IVW results:         ", nrow(ivw), "\n")
cat("MR summary:          ", MR_SUMMARY_FILE, "\n")
cat("Log file:            ", log_file, "\n")
