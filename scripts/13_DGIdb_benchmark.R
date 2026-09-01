#!/usr/bin/env Rscript

# ============================================================
# Framework_v1
# 13_DGIdb_benchmark.R
#
# External pharmacological validation of gene-drug modules
# against DGIdb.
#
# Benchmark unit:
#   Gene + NatMed drug category
#
# Accepted drug mappings:
#   REUSED_CPIC
#   MANUAL_ACCEPTED
#
# IMPORTANT
# ---------
# DGIdb is used here as broad external pharmacological support,
# not as a clinical PGx benchmark. No composite score is
# calculated. Interaction types, supporting sources and evidence
# are retained as descriptive annotations.
# ============================================================

options(stringsAsFactors = FALSE, warn = 1, scipen = 999)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
source(file.path(dirname(script_path), "..", "R", "bootstrap.R"))

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("The data.table package is required.", call. = FALSE)
}
library(data.table)

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

dgidb_file <- file.path(DIR_RAW, "DGIdb", "interactions.tsv")
mapping_file <- file.path(
  DIR_PROCESSED,
  "DGIdb_to_NatMed_DrugCategory_Mapping.tsv"
)
module_file <- file.path(
  PROJECT_ROOT,
  "results",
  "gene_drug_modules",
  "Gene_DrugModules.tsv"
)

output_dir <- file.path(PROJECT_ROOT, "results", "dgidb_benchmark")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

master_file <- file.path(
  output_dir,
  "DGIdbBenchmark_Master.tsv"
)

qc_file <- file.path(
  output_dir,
  "DGIdbBenchmark_QC.tsv"
)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

check_cols <- function(dt, cols, label) {
  miss <- setdiff(cols, names(dt))
  if (length(miss) > 0) {
    stop(
      label, " is missing required columns: ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )
  }
}

norm_drug <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("[[:space:]]+", " ", x)
  x
}

collapse_unique <- function(v, max_n = Inf) {
  v <- sort(unique(v[!is.na(v) & v != ""]))
  if (length(v) == 0) return(NA_character_)
  if (length(v) > max_n) {
    return(
      paste0(
        paste(v[seq_len(max_n)], collapse = ";"),
        ";...(+", length(v) - max_n, ")"
      )
    )
  }
  paste(v, collapse = ";")
}

safe_numeric_max <- function(x) {
  z <- suppressWarnings(as.numeric(x))
  if (all(is.na(z))) NA_real_ else max(z, na.rm = TRUE)
}

safe_numeric_median <- function(x) {
  z <- suppressWarnings(as.numeric(x))
  if (all(is.na(z))) NA_real_ else median(z, na.rm = TRUE)
}

# ------------------------------------------------------------
# Validate inputs
# ------------------------------------------------------------

required_files <- c(dgidb_file, mapping_file, module_file)
missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing required file(s):\n",
    paste(missing_files, collapse = "\n"),
    call. = FALSE
  )
}

cat("\n=== DGIdb external pharmacological benchmark ===\n")

# DGIdb snapshots may contain one or more leading metadata/comment
# lines, and the first data column is not guaranteed to be gene_name.
# Find the first non-empty, non-comment line and treat it as the header.
raw_lines <- readLines(dgidb_file, n = 100, warn = FALSE)

header_idx <- which(
  nzchar(trimws(raw_lines)) &
    !grepl("^\\s*#", raw_lines)
)[1]

if (is.na(header_idx)) {
  stop(
    "Could not locate a non-comment header line in DGIdb interactions.tsv.",
    call. = FALSE
  )
}

cat(
  "DGIdb header detected at file line: ",
  header_idx, "\n",
  sep = ""
)

dg <- fread(
  dgidb_file,
  sep = "\t",
  header = TRUE,
  skip = header_idx - 1L,
  na.strings = c("", "NA")
)

cat(
  "DGIdb columns detected: ",
  paste(names(dg), collapse = ", "),
  "\n"
)

mapping <- fread(
  mapping_file,
  sep = "\t",
  na.strings = c("", "NA")
)

modules <- fread(
  module_file,
  sep = "\t",
  na.strings = c("", "NA")
)

check_cols(
  dg,
  c(
    "gene_name",
    "drug_name",
    "drug_is_approved",
    "interaction_types",
    "interaction_source_db_name",
    "interaction_source_db_version",
    "interaction_score",
    "evidence_score"
  ),
  "DGIdb interactions.tsv"
)

check_cols(
  mapping,
  c(
    "DGIdb_Drug",
    "DGIdb_Drug_Normalized",
    "MappingStatus",
    "NatMed_DrugCategory",
    "NatMed_CategoryCode"
  ),
  "DGIdb_to_NatMed_DrugCategory_Mapping.tsv"
)

check_cols(
  modules,
  c(
    "AssignedGene",
    "DrugCategory",
    "CategoryCode",
    "UniqueSNPs",
    "UniqueMetabolites",
    "Tier1SNPs",
    "Tier1Metabolites",
    "Tier2SNPs",
    "Tier2Metabolites",
    "CodingSNPs",
    "SpliceRelatedSNPs",
    "HighImpactSNPs",
    "ModerateImpactSNPs",
    "MinMRFDR",
    "MaxPPH4",
    "MinM2P",
    "SNPs",
    "Metabolites",
    "CardiovascularATCGroup"
  ),
  "Gene_DrugModules.tsv"
)

# ------------------------------------------------------------
# Standardize raw DGIdb
# ------------------------------------------------------------

dg[, Gene := toupper(trimws(as.character(gene_name)))]
dg[, DGIdb_Drug_Normalized := norm_drug(drug_name)]

dg[, Approved :=
  tolower(trimws(as.character(drug_is_approved))) %in%
    c("true", "t", "1", "yes")
]

dg <- dg[
  !is.na(Gene) & Gene != "" &
    !is.na(DGIdb_Drug_Normalized) &
    DGIdb_Drug_Normalized != ""
]

approved <- dg[Approved == TRUE]

# ------------------------------------------------------------
# Accepted mappings only
# ------------------------------------------------------------

accepted_mapping <- mapping[
  MappingStatus %in% c("REUSED_CPIC", "MANUAL_ACCEPTED") &
    !is.na(NatMed_DrugCategory) &
    NatMed_DrugCategory != "" &
    !is.na(NatMed_CategoryCode) &
    NatMed_CategoryCode != ""
]

accepted_mapping[, DGIdb_Drug_Normalized :=
                   norm_drug(DGIdb_Drug_Normalized)]

if (nrow(accepted_mapping) == 0) {
  stop(
    "No accepted DGIdb mappings found. Expected MappingStatus ",
    "REUSED_CPIC or MANUAL_ACCEPTED.",
    call. = FALSE
  )
}

if (anyDuplicated(accepted_mapping$DGIdb_Drug_Normalized) > 0) {
  dup <- unique(
    accepted_mapping[
      duplicated(DGIdb_Drug_Normalized) |
        duplicated(DGIdb_Drug_Normalized, fromLast = TRUE),
      DGIdb_Drug_Normalized
    ]
  )

  stop(
    "Accepted DGIdb mapping contains duplicated normalized drug names: ",
    paste(head(dup, 20), collapse = ", "),
    call. = FALSE
  )
}

# ------------------------------------------------------------
# Individual DGIdb interaction rows -> accepted NatMed categories
# ------------------------------------------------------------

dg_mapped <- merge(
  approved,
  accepted_mapping[, .(
    DGIdb_Drug_Normalized,
    NatMed_DrugCategory,
    NatMed_CategoryCode,
    MappingStatus
  )],
  by = "DGIdb_Drug_Normalized",
  all = FALSE,
  sort = FALSE
)

if (nrow(dg_mapped) == 0) {
  stop(
    "No approved DGIdb interactions matched the accepted NatMed mapping.",
    call. = FALSE
  )
}

# ------------------------------------------------------------
# Collapse DGIdb evidence to Gene + NatMed drug category
# ------------------------------------------------------------

dg_category <- dg_mapped[
  ,
  .(
    DGIdb_InteractionRows = .N,
    DGIdb_UniqueDrugs = uniqueN(DGIdb_Drug_Normalized),
    DGIdb_Drugs = collapse_unique(drug_name, 50),
    DGIdb_InteractionTypes = collapse_unique(interaction_types, 50),
    DGIdb_UniqueInteractionTypes = uniqueN(
      interaction_types[!is.na(interaction_types) & interaction_types != ""]
    ),
    DGIdb_Sources = collapse_unique(interaction_source_db_name, 50),
    DGIdb_UniqueSources = uniqueN(
      interaction_source_db_name[
        !is.na(interaction_source_db_name) &
          interaction_source_db_name != ""
      ]
    ),
    DGIdb_SourceVersions = collapse_unique(
      interaction_source_db_version,
      50
    ),
    DGIdb_MaxEvidenceScore = safe_numeric_max(evidence_score),
    DGIdb_MedianEvidenceScore = safe_numeric_median(evidence_score),
    DGIdb_MaxInteractionScore = safe_numeric_max(interaction_score),
    DGIdb_MedianInteractionScore = safe_numeric_median(interaction_score),
    MappingStatuses = collapse_unique(MappingStatus, 10)
  ),
  by = .(
    Gene,
    NatMed_DrugCategory,
    NatMed_CategoryCode
  )
]

# One DGIdb benchmark unit per Gene + category.
if (anyDuplicated(
  dg_category[
    ,
    paste(
      Gene,
      NatMed_DrugCategory,
      NatMed_CategoryCode,
      sep = "||"
    )
  ]
) > 0) {
  stop(
    "Duplicate Gene + NatMed category rows remain in DGIdb collapse.",
    call. = FALSE
  )
}

# ------------------------------------------------------------
# Prepare our gene-drug modules
# ------------------------------------------------------------

module_key <- modules[, .(
  Gene = toupper(trimws(as.character(AssignedGene))),
  NatMed_DrugCategory = DrugCategory,
  NatMed_CategoryCode = CategoryCode,
  UniqueSNPs,
  UniqueMetabolites,
  Tier1SNPs,
  Tier1Metabolites,
  Tier2SNPs,
  Tier2Metabolites,
  CodingSNPs,
  SpliceRelatedSNPs,
  HighImpactSNPs,
  ModerateImpactSNPs,
  MinMRFDR,
  MaxPPH4,
  MinM2P,
  SNPs,
  Metabolites,
  CardiovascularATCGroup
)]

if (anyDuplicated(
  module_key[
    ,
    paste(
      Gene,
      NatMed_DrugCategory,
      NatMed_CategoryCode,
      sep = "||"
    )
  ]
) > 0) {
  stop(
    "Gene_DrugModules.tsv contains duplicated gene-category modules.",
    call. = FALSE
  )
}

# ------------------------------------------------------------
# Benchmark universe = mapped DGIdb gene-category units
# ------------------------------------------------------------
# Left join preserves every DGIdb-supported benchmark unit.
# Recovered indicates whether the same Gene + NatMed category
# exists in our metabolically anchored gene-drug modules.

master <- merge(
  dg_category,
  module_key,
  by = c(
    "Gene",
    "NatMed_DrugCategory",
    "NatMed_CategoryCode"
  ),
  all.x = TRUE,
  sort = FALSE
)

master[, Recovered := !is.na(UniqueSNPs)]

# Descriptive mechanism flags based on DGIdb interaction-type text.
itype <- tolower(
  fifelse(
    is.na(master$DGIdb_InteractionTypes),
    "",
    master$DGIdb_InteractionTypes
  )
)

master[, DGIdb_InhibitoryEvidence := grepl(
  "inhibitor|inhibitory|antagonist|blocker|negative modulator",
  itype
)]

master[, DGIdb_ActivatingEvidence := grepl(
  "activator|agonist|positive modulator|stimulator",
  itype
)]

master[, DGIdb_TransportEvidence := grepl(
  "transporter|transport|substrate",
  itype
)]

master[, DGIdb_TargetLikeEvidence := grepl(
  "inhibitor|antagonist|agonist|activator|modulator|blocker|binder|ligand",
  itype
)]

setorder(
  master,
  -Recovered,
  -DGIdb_UniqueSources,
  -DGIdb_UniqueDrugs,
  Gene,
  NatMed_DrugCategory
)

fwrite(
  master,
  master_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# ------------------------------------------------------------
# QC / summary
# ------------------------------------------------------------

safe_pct <- function(n, d) {
  if (d == 0) return(NA_real_)
  100 * n / d
}

n_benchmark <- nrow(master)
n_recovered <- master[Recovered == TRUE, .N]


qc <- data.table(
  Metric = c(
    "Raw DGIdb interaction rows",
    "Approved DGIdb interaction rows",
    "Accepted mapped DGIdb drugs",
    "Mapped approved DGIdb interaction rows",
    "DGIdb Gene-NatMed category benchmark units",
    "DGIdb benchmark units recovered",
    "DGIdb category-level recovery percent",
    "Unique DGIdb genes in benchmark",
    "Unique NatMed categories in DGIdb benchmark",
    "Recovered unique genes",
    "Recovered NatMed categories",
    "Recovered units with >=2 DGIdb drugs",
    "Recovered units with >=2 DGIdb sources",
    "Recovered units with coding SNP support",
    "Recovered units with splice-related SNP support",
    "Recovered units with HIGH-impact SNP support",
    "Recovered units with DGIdb inhibitory evidence",
    "Recovered units with DGIdb activating evidence",
    "Recovered units with DGIdb transport evidence",
    "Recovered units with DGIdb target-like evidence"
  ),
  Value = c(
    as.character(nrow(dg)),
    as.character(nrow(approved)),
    as.character(uniqueN(accepted_mapping$DGIdb_Drug_Normalized)),
    as.character(nrow(dg_mapped)),
    as.character(n_benchmark),
    as.character(n_recovered),
    sprintf("%.2f", safe_pct(n_recovered, n_benchmark)),
    as.character(uniqueN(master$Gene)),
    as.character(uniqueN(master$NatMed_DrugCategory)),
    as.character(uniqueN(master[Recovered == TRUE, Gene])),
    as.character(uniqueN(
      master[Recovered == TRUE, NatMed_DrugCategory]
    )),
    as.character(
      master[
        Recovered == TRUE &
          DGIdb_UniqueDrugs >= 2,
        .N
      ]
    ),
    as.character(
      master[
        Recovered == TRUE &
          DGIdb_UniqueSources >= 2,
        .N
      ]
    ),
    as.character(
      master[
        Recovered == TRUE &
          CodingSNPs > 0,
        .N
      ]
    ),
    as.character(
      master[
        Recovered == TRUE &
          SpliceRelatedSNPs > 0,
        .N
      ]
    ),
    as.character(
      master[
        Recovered == TRUE &
          HighImpactSNPs > 0,
        .N
      ]
    ),
    as.character(
      master[
        Recovered == TRUE &
          DGIdb_InhibitoryEvidence == TRUE,
        .N
      ]
    ),
    as.character(
      master[
        Recovered == TRUE &
          DGIdb_ActivatingEvidence == TRUE,
        .N
      ]
    ),
    as.character(
      master[
        Recovered == TRUE &
          DGIdb_TransportEvidence == TRUE,
        .N
      ]
    ),
    as.character(
      master[
        Recovered == TRUE &
          DGIdb_TargetLikeEvidence == TRUE,
        .N
      ]
    )
  )
)

fwrite(
  qc,
  qc_file,
  sep = "\t",
  quote = FALSE
)

cat("\n============================================================\n")
cat("DGIdb benchmark completed successfully\n")
cat("============================================================\n\n")
print(qc)

cat("\nSelected recovered modules:\n")

focus_genes <- c(
  "HMGCR", "SLCO1B1", "ADRB1", "ABCB1",
  "CETP", "LDLR", "LPA", "APOE",
  "CYP2C19", "CYP2D6", "CYP2C9"
)

print(
  master[
    Recovered == TRUE &
      Gene %in% focus_genes,
    .(
      Gene,
      NatMed_DrugCategory,
      DGIdb_UniqueDrugs,
      DGIdb_UniqueSources,
      DGIdb_InteractionTypes,
      UniqueSNPs,
      UniqueMetabolites,
      CodingSNPs,
      SpliceRelatedSNPs,
      DGIdb_InhibitoryEvidence,
      DGIdb_TransportEvidence
    )
  ]
)

cat(
  "\nOutputs:\n",
  "  ", master_file, "\n",
  "  ", qc_file, "\n\n",
  sep = ""
)

cat(
  "INTERPRETATION NOTE:\n",
  "DGIdb support is broad pharmacological interaction evidence and\n",
  "must not be interpreted as clinical PGx actionability. Recovery\n",
  "is evaluated at Gene + NatMed drug-category resolution.\n"
)
