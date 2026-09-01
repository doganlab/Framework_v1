# Framework_v1 shared bootstrap --------------------------------------------
get_framework_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (!length(hit)) return(NULL)
  normalizePath(sub("^--file=", "", hit[[1]]), winslash = "/", mustWork = TRUE)
}
resolve_framework_root <- function() {
  script <- get_framework_script_path()
  candidate <- if (is.null(script)) getwd() else file.path(dirname(script), "..")
  normalizePath(candidate, winslash = "/", mustWork = TRUE)
}
PROJECT_ROOT <- resolve_framework_root()
setwd(PROJECT_ROOT)
source(file.path(PROJECT_ROOT, "config", "framework_config.R"))
