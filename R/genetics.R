# ============================================================
# Framework_v1
# Genetics helper functions
# ============================================================

#--------------------------------------------------------------
# Create standardized variant identifier
#--------------------------------------------------------------

create_variant_id <- function(chr, bp, effect_allele = NA, other_allele = NA){
  
  if(all(is.na(effect_allele)) || all(is.na(other_allele))){
    
    paste0("chr", chr, ":", bp)
    
  } else {
    
    paste0(
      "chr", chr, ":", bp, ":",
      effect_allele, ":",
      other_allele
    )
    
  }
  
}

#--------------------------------------------------------------
# Create locus ID
#--------------------------------------------------------------

create_locus_id <- function(chr, bp, window_kb = 250){
  
  locus_start <- floor(as.numeric(bp) / (window_kb * 1000)) * window_kb * 1000
  locus_end <- locus_start + window_kb * 1000 - 1
  
  paste0(
    "chr", chr, ":",
    locus_start, "-",
    locus_end
  )
  
}

#--------------------------------------------------------------
# Infer direction of effect
#--------------------------------------------------------------

infer_direction <- function(mr_beta){
  
  direction <- rep(NA_character_, length(mr_beta))
  
  direction[mr_beta > 0] <- "Risk_increasing"
  direction[mr_beta < 0] <- "Risk_decreasing"
  direction[mr_beta == 0] <- "Neutral"
  
  direction
  
}

#--------------------------------------------------------------
# Standardize chromosome names
#--------------------------------------------------------------

standardize_chr <- function(chr){
  
  chr <- as.character(chr)
  chr <- gsub("^chr", "", chr, ignore.case = TRUE)
  chr
  
}

#--------------------------------------------------------------
# Standardize SNP column
#--------------------------------------------------------------

standardize_snp <- function(snp){
  
  as.character(snp)
  
}