# ============================================================
# Framework_v1
# MR clumping functions
# ============================================================

library(data.table)

#--------------------------------------------------------------
# Run PLINK LD clumping for one exposure file
#--------------------------------------------------------------

run_plink_clump <- function(
    input_file,
    output_prefix,
    plink_bin,
    ld_ref_prefix,
    clump_p1 = 5e-8,
    clump_p2 = 5e-8,
    clump_r2 = 0.001,
    clump_kb = 10000
){
  
  check_file_exists(input_file)
  
  cmd <- paste(
    shQuote(plink_bin),
    "--bfile", shQuote(ld_ref_prefix),
    "--clump", shQuote(input_file),
    "--clump-snp-field SNP",
    "--clump-field P",
    "--clump-p1", clump_p1,
    "--clump-p2", clump_p2,
    "--clump-r2", clump_r2,
    "--clump-kb", clump_kb,
    "--out", shQuote(output_prefix)
  )
  
  log_message(paste0("Running PLINK clumping: ", basename(input_file)))
  
  status <- system(cmd)
  
  if(status != 0){
    stop_framework(
      paste0("PLINK clumping failed for: ", input_file)
    )
  }
  
  clumped_file <- paste0(output_prefix, ".clumped")
  
  check_file_exists(clumped_file)
  
  clumped_file
  
}

#--------------------------------------------------------------
# Run PLINK LD clumping for all prepared exposure files
#--------------------------------------------------------------

run_plink_clump_batch <- function(
    input_dir,
    output_dir,
    plink_bin,
    ld_ref_prefix,
    pattern = "_GWS_for_clumping.tsv$",
    clump_p1 = 5e-8,
    clump_p2 = 5e-8,
    clump_r2 = 0.001,
    clump_kb = 10000
){
  
  check_dir_exists(input_dir)
  create_directory(output_dir)
  
  files <- list.files(
    input_dir,
    pattern = pattern,
    full.names = TRUE
  )
  
  if(length(files) == 0){
    stop_framework(
      paste0("No clumping input files found in: ", input_dir)
    )
  }
  
  result_list <- list()
  
  for(file in files){
    
    metabolite <- sub(pattern, "", basename(file))
    
    output_prefix <- file.path(output_dir, metabolite)
    
    clumped_file <- run_plink_clump(
      input_file = file,
      output_prefix = output_prefix,
      plink_bin = plink_bin,
      ld_ref_prefix = ld_ref_prefix,
      clump_p1 = clump_p1,
      clump_p2 = clump_p2,
      clump_r2 = clump_r2,
      clump_kb = clump_kb
    )
    
    result_list[[metabolite]] <- data.table(
      Metabolite = metabolite,
      InputFile = file,
      ClumpedFile = clumped_file
    )
    
  }
  
  rbindlist(result_list, fill = TRUE)
  
}

#--------------------------------------------------------------
# Summarize clumped instruments
#--------------------------------------------------------------

summarize_clumped_instruments <- function(clump_dir){
  
  check_dir_exists(clump_dir)
  
  files <- list.files(
    clump_dir,
    pattern = "\\.clumped$",
    full.names = TRUE
  )
  
  if(length(files) == 0){
    stop_framework(
      paste0("No .clumped files found in: ", clump_dir)
    )
  }
  
  res <- rbindlist(
    lapply(files, function(f){
      
      met <- sub("\\.clumped$", "", basename(f))
      x <- fread(f, fill = TRUE)
      
      data.table(
        Metabolite = met,
        n_clumped_snps = nrow(x)
      )
      
    }),
    fill = TRUE
  )
  
  setorder(res, -n_clumped_snps)
  
  res
  
}