# ============================================================
# Framework_v1
# MR input preparation functions -- streaming mGWAS edition
# ============================================================

library(data.table)

load_candidate_metabolites <- function(candidate_file){
  check_file_exists(candidate_file)
  candidates <- fread(candidate_file)
  required_input_cols <- c("metabolite", "study_accession", "PP.H4")
  check_columns(candidates, required_input_cols)
  setnames(candidates,
           old = c("metabolite", "study_accession", "PP.H4"),
           new = c("Metabolite", "StudyAccession", "PP_H4"))
  candidates[, Metabolite := trimws(as.character(Metabolite))]
  candidates[, StudyAccession := trimws(as.character(StudyAccession))]
  candidates[, PP_H4 := as.numeric(PP_H4)]
  candidates <- candidates[
    !is.na(Metabolite) & Metabolite != "" &
      !is.na(StudyAccession) & StudyAccession != ""
  ]
  candidates <- unique(candidates, by = c("Metabolite", "StudyAccession"))
  candidates
}

locate_candidate_mgwas_files <- function(candidates, mgwas_dir){
  check_dir_exists(mgwas_dir)
  check_columns(candidates, c("Metabolite", "StudyAccession"))
  files <- list.files(mgwas_dir, pattern = "\\.tsv\\.gz$", full.names = TRUE)
  candidates[, mGWAS_file := NA_character_]
  for(i in seq_len(nrow(candidates))){
    accession <- candidates$StudyAccession[i]
    hit <- files[grepl(accession, basename(files), fixed = TRUE)]
    if(length(hit) == 1){
      candidates$mGWAS_file[i] <- hit
    } else if(length(hit) > 1){
      warn_framework(paste0("Multiple mGWAS files found for accession ", accession,
                            ". Using the first match."))
      candidates$mGWAS_file[i] <- hit[1]
    }
  }
  candidates
}

# -----------------------------------------------------------------------------
# Stream only genome-wide significant records from a large gzip-compressed
# GWAS Catalog file. The complete ~3 GB table is never materialized in R RAM.
# -----------------------------------------------------------------------------
extract_gws_instruments_streaming <- function(mgwas_file, p_threshold = 5e-8){
  check_file_exists(mgwas_file)
  if(!is.numeric(p_threshold) || length(p_threshold) != 1L ||
     is.na(p_threshold) || p_threshold <= 0 || p_threshold >= 1){
    stop_framework("p_threshold must be a single number between 0 and 1.")
  }

  neglog_threshold <- -log10(p_threshold)
  f <- shQuote(normalizePath(mgwas_file, mustWork = TRUE))

  awk_program <- paste0(
    "'BEGIN{FS=\"\\t\";OFS=\"\\t\"} ",
    "NR==1{",
    "for(i=1;i<=NF;i++)h[$i]=i;",
    "req[1]=\"rsid\";req[2]=\"chromosome\";req[3]=\"base_pair_location\";",
    "req[4]=\"effect_allele\";req[5]=\"other_allele\";req[6]=\"beta\";",
    "req[7]=\"standard_error\";req[8]=\"effect_allele_frequency\";",
    "req[9]=\"neg_log_10_p_value\";",
    "for(j=1;j<=9;j++)if(!(req[j] in h)){print \"Missing column: \" req[j] > \"/dev/stderr\";exit 42};",
    "print \"SNP\",\"CHR\",\"BP\",\"effect_allele\",\"other_allele\",",
    "\"beta\",\"se\",\"eaf\",\"neg_log10_p\";next}",
    " $(h[\"neg_log_10_p_value\"])+0 > ", format(neglog_threshold, digits = 16),
    " {print $(h[\"rsid\"]),$(h[\"chromosome\"]),$(h[\"base_pair_location\"]),",
    "$(h[\"effect_allele\"]),$(h[\"other_allele\"]),$(h[\"beta\"]),",
    "$(h[\"standard_error\"]),$(h[\"effect_allele_frequency\"]),",
    "$(h[\"neg_log_10_p_value\"])}'"
  )

  cmd <- paste("gzip -cd", f, "| awk", awk_program)
  x <- tryCatch(
    fread(cmd = cmd, sep = "\t", header = TRUE, showProgress = FALSE),
    error = function(e){
      stop_framework(paste0("Streaming GWS extraction failed for ", basename(mgwas_file),
                            ": ", conditionMessage(e)))
    }
  )

  if(nrow(x) == 0L) return(x)

  x[, effect_allele := toupper(trimws(as.character(effect_allele)))]
  x[, other_allele := toupper(trimws(as.character(other_allele)))]
  x[, `:=`(
    SNP = trimws(as.character(SNP)),
    CHR = suppressWarnings(as.integer(CHR)),
    BP = suppressWarnings(as.integer(BP)),
    beta = suppressWarnings(as.numeric(beta)),
    se = suppressWarnings(as.numeric(se)),
    eaf = suppressWarnings(as.numeric(eaf)),
    neg_log10_p = suppressWarnings(as.numeric(neg_log10_p))
  )]
  x[, P := 10^(-neg_log10_p)]

  x <- x[
    !is.na(SNP) & SNP != "" & !is.na(CHR) & !is.na(BP) &
      !is.na(effect_allele) & effect_allele != "" &
      !is.na(other_allele) & other_allele != "" &
      !is.na(beta) & !is.na(se) & se > 0 & !is.na(P) & P < p_threshold
  ]
  setorder(x, P)
  unique(x, by = "SNP")
}

# Backward-compatible name. For these large Nightingale files this now streams
# only GWS rows instead of reading all ~96 million variants into R.
read_mgwas_standard <- function(mgwas_file, p_threshold = 5e-8){
  extract_gws_instruments_streaming(mgwas_file, p_threshold = p_threshold)
}

select_gws_instruments <- function(mgwas, p_threshold = 5e-8){
  required_cols <- c("SNP", "CHR", "BP", "effect_allele", "other_allele",
                     "beta", "se", "eaf", "P")
  check_columns(mgwas, required_cols)
  mgwas[P < p_threshold]
}

prepare_clumping_input <- function(instruments){
  check_columns(instruments, c("SNP", "P"))
  instruments[, .(SNP, P)]
}

prepare_all_clumping_inputs <- function(candidates, output_dir, p_threshold = 5e-8){
  create_directory(output_dir)
  summary_list <- list()

  for(i in seq_len(nrow(candidates))){
    met <- candidates$Metabolite[i]
    mgwas_file <- candidates$mGWAS_file[i]
    log_message(paste0("Preparing clumping input for: ", met))

    if(is.na(mgwas_file) || mgwas_file == ""){
      warn_framework(paste0("No mGWAS file found for: ", met))
      summary_list[[met]] <- data.table(
        Metabolite = met, mGWAS_file = NA_character_, n_gws_snps = NA_integer_,
        gws_instrument_file = NA_character_, clumping_input_file = NA_character_,
        status = "mGWAS_not_found"
      )
      next
    }

    safe_met <- gsub("[^A-Za-z0-9_.-]+", "_", met)
    gws_file <- file.path(output_dir, paste0(safe_met, "_GWS_instruments.tsv.gz"))
    clump_file <- file.path(output_dir, paste0(safe_met, "_GWS_for_clumping.tsv"))

    # Resume-safe: reuse a non-empty GWS cache generated previously.
    if(file.exists(gws_file) && file.info(gws_file)$size > 20){
      log_message(paste0("Reusing cached GWS instruments: ", basename(gws_file)))
      instruments <- fread(gws_file, showProgress = FALSE)
      # Guard against a cache produced with a different/looser threshold.
      if("P" %in% names(instruments)) instruments <- instruments[P < p_threshold]
    } else {
      instruments <- extract_gws_instruments_streaming(
        mgwas_file = mgwas_file,
        p_threshold = p_threshold
      )
      if(nrow(instruments) > 0L){
        fwrite(instruments, gws_file, sep = "\t", compress = "gzip")
      }
    }

    if(nrow(instruments) == 0L){
      warn_framework(paste0("No GWS instruments for: ", met))
      summary_list[[met]] <- data.table(
        Metabolite = met, mGWAS_file = mgwas_file, n_gws_snps = 0L,
        gws_instrument_file = NA_character_, clumping_input_file = NA_character_,
        status = "no_gws_instruments"
      )
      next
    }

    clump_input <- prepare_clumping_input(instruments)
    fwrite(clump_input, clump_file, sep = "\t")

    summary_list[[met]] <- data.table(
      Metabolite = met,
      mGWAS_file = mgwas_file,
      n_gws_snps = nrow(instruments),
      gws_instrument_file = gws_file,
      clumping_input_file = clump_file,
      status = "ok"
    )

    rm(instruments, clump_input)
    invisible(gc(FALSE))
  }

  rbindlist(summary_list, fill = TRUE)
}
