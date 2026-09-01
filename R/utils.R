# ============================================================
# Framework_v1
# Utility functions
# ============================================================

#--------------------------------------------------------------
# Log message
#--------------------------------------------------------------

log_message <- function(message){
  
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  
  cat(
    paste0("[", timestamp, "] ", message, "\n")
  )
  
}

#--------------------------------------------------------------
# Stop with formatted message
#--------------------------------------------------------------

stop_framework <- function(message){
  
  stop(
    paste0("[Framework_v1 ERROR] ", message),
    call. = FALSE
  )
  
}

#--------------------------------------------------------------
# Warn with formatted message
#--------------------------------------------------------------

warn_framework <- function(message){
  
  warning(
    paste0("[Framework_v1 WARNING] ", message),
    call. = FALSE
  )
  
}

#--------------------------------------------------------------
# Check file exists
#--------------------------------------------------------------

check_file_exists <- function(file){
  
  if(!file.exists(file)){
    stop_framework(
      paste0("File not found: ", file)
    )
  }
  
  invisible(TRUE)
  
}

#--------------------------------------------------------------
# Check directory exists
#--------------------------------------------------------------

check_dir_exists <- function(path){
  
  if(!dir.exists(path)){
    stop_framework(
      paste0("Directory not found: ", path)
    )
  }
  
  invisible(TRUE)
  
}