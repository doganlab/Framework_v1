# ============================================================
# Framework_v1
# IO functions
# ============================================================

library(data.table)

#--------------------------------------------------------------
# Read TSV
#--------------------------------------------------------------

read_tsv_safe <- function(file){
  
  if(!file.exists(file)){
    stop("File not found: ", file)
  }
  
  fread(
    file,
    sep="\t",
    header=TRUE,
    data.table=FALSE
  )
  
}

#--------------------------------------------------------------
# Write TSV
#--------------------------------------------------------------

write_tsv_safe <- function(df, file){
  
  fwrite(
    df,
    file,
    sep="\t",
    quote=FALSE,
    row.names=FALSE,
    na="NA"
  )
  
}

#--------------------------------------------------------------
# Check required columns
#--------------------------------------------------------------

check_columns <- function(df, required_columns){
  
  missing <- setdiff(required_columns, colnames(df))
  
  if(length(missing)>0){
    
    stop(
      "Missing columns:\n",
      paste(missing, collapse="\n")
    )
    
  }
  
  invisible(TRUE)
  
}


create_directory <- function(path){
  
  if(!dir.exists(path)){
    
    dir.create(
      path,
      recursive = TRUE,
      showWarnings = FALSE
    )
    
  }
  
}

