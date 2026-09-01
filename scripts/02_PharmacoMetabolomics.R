#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, scipen = 999)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)

source(file.path(dirname(script_path), "..", "R", "bootstrap.R"))
source(file.path(PROJECT_ROOT, "R", "io.R"))
source(file.path(PROJECT_ROOT, "R", "utils.R"))
source(file.path(PROJECT_ROOT, "R", "pharmacometabolomics.R"))

required_files <- c(
  MR_SUMMARY_FILE,
  CANDIDATE_METABOLITE_FILE,
  NATMED_M2_ASSOCIATION_FILE,
  NATMED_DRUG_CATEGORY_FILE
)

missing <- required_files[!file.exists(required_files)]

if (length(missing) > 0) {
  stop(
    "Missing required input file(s):\n",
    paste(missing, collapse = "\n"),
    call. = FALSE
  )
}

create_directory(DIR_PROCESSED)
create_directory(DIR_LOGS)

cat("\nBuilding pharmacometabolomic evidence table...\n")

mr <- data.table::fread(MR_SUMMARY_FILE)
mapping <- data.table::fread(CANDIDATE_METABOLITE_FILE)
assoc <- data.table::fread(NATMED_M2_ASSOCIATION_FILE)
drugs <- data.table::fread(NATMED_DRUG_CATEGORY_FILE)

cat("MR rows:                  ", nrow(mr), "\n", sep = "")
cat("Mapped candidate traits:  ", data.table::uniqueN(mapping$metabolite), "\n", sep = "")
cat("NatMed association rows:  ", nrow(assoc), "\n", sep = "")
cat("NatMed drug categories:   ", data.table::uniqueN(drugs$DrugCategory), "\n", sep = "")

evidence <- build_pharmacometabolomic_evidence(
  mr_summary = mr,
  candidate_mapping = mapping,
  m2_associations = assoc,
  drug_categories = drugs,
  fallback_threshold = NATMED_M2_FALLBACK_P_THRESHOLD
)

if (nrow(evidence) == 0) {
  stop(
    "No MR metabolite could be joined to Nature Medicine Suppl_T4 associations.",
    call. = FALSE
  )
}

# For the currently validated source files:
# 102 matched traits x 87 drug categories = 8,874 rows.
expected_rows <- data.table::uniqueN(evidence$Metabolite) *
  data.table::uniqueN(drugs$DrugCategory)

if (nrow(evidence) != expected_rows) {
  stop(
    paste0(
      "Unexpected pharmacometabolomic row count. Observed: ",
      nrow(evidence),
      "; expected from represented metabolites x drug categories: ",
      expected_rows,
      "."
    ),
    call. = FALSE
  )
}

write_tsv_safe(
  evidence,
  PHARMACOMETABOLOMICS_FILE
)

summary <- data.table::data.table(
  Metric = c(
    "Rows",
    "MR metabolites represented",
    "Drug categories represented",
    "M2 significant associations",
    "M2 significant drug categories",
    "M2 significant and risk-opposing associations",
    "M2 significant, risk-opposing, MR-FDR significant associations",
    "Tier1 robust associations",
    "Tier2 high-confidence associations",
    "Directionally discordant associations",
    "Cardiovascular ATC rows among primary associations"
  ),
  Value = c(
    nrow(evidence),
    data.table::uniqueN(evidence$Metabolite),
    data.table::uniqueN(evidence$DrugCategory),
    evidence[M2Significant %in% TRUE, .N],
    data.table::uniqueN(
      evidence[M2Significant %in% TRUE, DrugCategory]
    ),
    evidence[
      M2Significant %in% TRUE &
        RiskOpposingDirection %in% TRUE,
      .N
    ],
    evidence[
      PrimaryPharmacoGeneticEvidence %in% TRUE,
      .N
    ],
    evidence[EvidenceTier == "Tier1_Robust", .N],
    evidence[EvidenceTier == "Tier2_HighConfidence", .N],
    evidence[EvidenceTier == "DirectionallyDiscordant", .N],
    evidence[
      PrimaryPharmacoGeneticEvidence %in% TRUE &
        CardiovascularATCGroup %in% TRUE,
      .N
    ]
  )
)

write_tsv_safe(
  summary,
  PHARMACOMETABOLOMICS_SUMMARY_FILE
)

cat("\nPharmacometabolomic integration completed.\n")
print(summary)

cat(
  "\nOutput: ",
  PHARMACOMETABOLOMICS_FILE,
  "\nSummary: ",
  PHARMACOMETABOLOMICS_SUMMARY_FILE,
  "\n\n",
  sep = ""
)
