#!/usr/bin/env Rscript

# ============================================================
# Framework_v1
# 12_DGIdb_prepare_mapping.R
#
# Stage 1 of DGIdb external pharmacological validation.
#
# Purpose:
#   1. Read the raw DGIdb interactions.tsv snapshot.
#   2. Restrict the mapping universe to approved drugs.
#   3. Reuse the manually reviewed CPIC -> NatMed mapping wherever
#      an identical normalized drug name is available.
#   4. Identify additional DGIdb approved drugs that are plausible
#      candidates for mapping to one of the NatMed drug categories.
#   5. Write a conservative mapping-review table.
#
# IMPORTANT:
#   This script DOES NOT perform the gene-drug benchmark yet.
#   It prepares the drug-category mapping that must be reviewed first.
# ============================================================

options(stringsAsFactors = FALSE, warn = 1, scipen = 999)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
source(file.path(dirname(script_path), "..", "R", "bootstrap.R"))

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("The data.table package is required.", call. = FALSE)
}
library(data.table)

dgidb_file <- file.path(DIR_RAW, "DGIdb", "interactions.tsv")
cpic_mapping_file <- file.path(
  DIR_PROCESSED,
  "CPIC_to_NatMed_DrugCategory_Mapping.tsv"
)
module_file <- file.path(
  PROJECT_ROOT,
  "results",
  "gene_drug_modules",
  "Gene_DrugModules.tsv"
)

out_file <- file.path(
  DIR_PROCESSED,
  "DGIdb_to_NatMed_DrugCategory_Mapping.tsv"
)
qc_file <- file.path(
  DIR_PROCESSED,
  "DGIdb_MappingPreparation_QC.tsv"
)

required <- c(dgidb_file, cpic_mapping_file, module_file)
missing <- required[!file.exists(required)]
if (length(missing)) {
  stop(
    "Missing required file(s):\n",
    paste(missing, collapse = "\n"),
    call. = FALSE
  )
}

cat("\n=== DGIdb -> NatMed mapping preparation ===\n")

# fread automatically ignores the two leading DGIdb comment lines.
dg <- fread(
  dgidb_file,
  sep = "\t",
  header = TRUE,
  na.strings = c("", "NA"),
  comment.char = "#"
)

cpic_map <- fread(
  cpic_mapping_file,
  sep = "\t",
  na.strings = c("", "NA")
)

modules <- fread(
  module_file,
  sep = "\t",
  na.strings = c("", "NA")
)

check_cols <- function(dt, cols, label) {
  miss <- setdiff(cols, names(dt))
  if (length(miss)) {
    stop(
      label, " is missing required columns: ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )
  }
}

check_cols(
  dg,
  c(
    "gene_name",
    "drug_name",
    "drug_is_approved",
    "interaction_source_db_name",
    "interaction_source_db_version",
    "interaction_types",
    "interaction_score",
    "evidence_score"
  ),
  "DGIdb interactions.tsv"
)

check_cols(
  cpic_map,
  c(
    "CPIC_Drug",
    "CPIC_Drug_Normalized",
    "NatMed_DrugCategory",
    "NatMed_CategoryCode",
    "MappingStatus",
    "MappingBasis",
    "MappingNotes"
  ),
  "CPIC_to_NatMed_DrugCategory_Mapping.tsv"
)

check_cols(
  modules,
  c("DrugCategory", "CategoryCode"),
  "Gene_DrugModules.tsv"
)

norm_drug <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("[[:space:]]+", " ", x)
  x
}

# Standardize approved flag robustly.
dg[, Approved := tolower(trimws(as.character(drug_is_approved))) %in%
     c("true", "t", "1", "yes")]

dg[, DrugNormalized := norm_drug(drug_name)]
dg[, gene_name := toupper(trimws(as.character(gene_name)))]

dg <- dg[
  !is.na(drug_name) & drug_name != "" &
  !is.na(gene_name) & gene_name != ""
]

approved <- dg[Approved == TRUE]

# Build one row per approved DGIdb drug with useful review context.
drug_summary <- approved[, .(
  DGIdb_InteractionRows = .N,
  DGIdb_UniqueGenes = uniqueN(gene_name),
  DGIdb_InteractionTypes = paste(
    sort(unique(na.omit(interaction_types))),
    collapse = ";"
  ),
  DGIdb_Sources = paste(
    sort(unique(na.omit(interaction_source_db_name))),
    collapse = ";"
  ),
  DGIdb_MaxEvidenceScore = {
    z <- suppressWarnings(as.numeric(evidence_score))
    if (all(is.na(z))) NA_real_ else max(z, na.rm = TRUE)
  },
  DGIdb_MaxInteractionScore = {
    z <- suppressWarnings(as.numeric(interaction_score))
    if (all(is.na(z))) NA_real_ else max(z, na.rm = TRUE)
  }
), by = .(
  DGIdb_Drug = drug_name,
  DGIdb_Drug_Normalized = DrugNormalized
)]

# CPIC mapping: only EXACT / CLASS are automatically reusable.
cpic_reuse <- copy(cpic_map)
cpic_reuse[, CPIC_Drug_Normalized :=
             norm_drug(CPIC_Drug_Normalized)]

cpic_reuse <- cpic_reuse[
  MappingStatus %in% c("EXACT", "CLASS") &
  !is.na(NatMed_DrugCategory) &
  NatMed_DrugCategory != "",
  .(
    DGIdb_Drug_Normalized = CPIC_Drug_Normalized,
    CPIC_SourceDrug = CPIC_Drug,
    NatMed_DrugCategory,
    NatMed_CategoryCode,
    CPIC_MappingStatus = MappingStatus,
    CPIC_MappingBasis = MappingBasis,
    CPIC_MappingNotes = MappingNotes
  )
]

if (anyDuplicated(cpic_reuse$DGIdb_Drug_Normalized)) {
  stop(
    "Accepted CPIC mapping has duplicated normalized drug names.",
    call. = FALSE
  )
}

mapping <- merge(
  drug_summary,
  cpic_reuse,
  by = "DGIdb_Drug_Normalized",
  all.x = TRUE,
  sort = FALSE
)

mapping[, MappingStatus := fifelse(
  !is.na(NatMed_DrugCategory),
  "REUSED_CPIC",
  "UNREVIEWED"
)]

mapping[, MappingBasis := fifelse(
  MappingStatus == "REUSED_CPIC",
  paste0(
    "Reused accepted CPIC mapping (",
    CPIC_MappingStatus,
    ")"
  ),
  "No accepted CPIC mapping available"
)]

mapping[, MappingNotes := fifelse(
  MappingStatus == "REUSED_CPIC",
  fifelse(
    is.na(CPIC_MappingNotes),
    "",
    CPIC_MappingNotes
  ),
  ""
)]

# ----------------------------------------------------------------
# Conservative candidate flagging for manual review.
#
# We do NOT automatically assign new categories here. Instead, we
# flag unreviewed DGIdb drugs whose names plausibly indicate one of
# the 29 NatMed categories. This keeps mapping auditable.
# ----------------------------------------------------------------

category_rules <- data.table(
  Pattern = c(
    "statin$|atorvastatin|simvastatin|lovastatin|fluvastatin|pravastatin|rosuvastatin|pitavastatin|cerivastatin",
    "metformin",
    "clopidogrel",
    "aspirin|acetylsalicylic",
    "omeprazole|esomeprazole|lansoprazole|dexlansoprazole|pantoprazole|rabeprazole",
    "captopril|enalapril|lisinopril|ramipril|perindopril|quinapril|fosinopril|benazepril|trandolapril",
    "atenolol|metoprolol|bisoprolol|betaxolol|acebutolol|esmolol|nebivolol",
    "furosemide|bumetanide|torsemide",
    "hydrochlorothiazide|chlorthalidone|indapamide",
    "spironolactone|eplerenone",
    "salbutamol|albuterol|salmeterol|formoterol",
    "gliclazide|glimepiride|glipizide|glyburide|glibenclamide|chlorpropamide|tolazamide|tolbutamide"
  ),
  CandidateHint = c(
    "Statin category - review lipophilic vs hydrophilic",
    "Metformin",
    "Clopidogrel",
    "Acetylsalicylic acid",
    "Proton pump inhibitors",
    "ACE inhibitors - plain",
    "Selective beta blockers",
    "High-ceiling diuretics",
    "Low-ceiling diuretics",
    "Potassium-sparing agents",
    "Selective beta-2-adrenoreceptor agonists",
    "Sulfonamides-urea derivatives"
  )
)

mapping[, CandidateNatMedHint := NA_character_]

for (i in seq_len(nrow(category_rules))) {
  hit <- mapping[
    MappingStatus == "UNREVIEWED" &
    grepl(
      category_rules$Pattern[i],
      DGIdb_Drug_Normalized,
      perl = TRUE
    ),
    which = TRUE
  ]

  if (length(hit)) {
    mapping[hit, CandidateNatMedHint :=
              fifelse(
                is.na(CandidateNatMedHint),
                category_rules$CandidateHint[i],
                paste(
                  CandidateNatMedHint,
                  category_rules$CandidateHint[i],
                  sep = ";"
                )
              )]
  }
}

mapping[, ReviewPriority := fifelse(
  MappingStatus == "REUSED_CPIC",
  "ACCEPTED",
  fifelse(
    !is.na(CandidateNatMedHint),
    "REVIEW",
    "LOW"
  )
)]

# Ensure reused category/code combinations exist in the actual modules.
known <- unique(modules[, .(
  NatMed_DrugCategory = DrugCategory,
  NatMed_CategoryCode = CategoryCode
)])

bad <- mapping[
  MappingStatus == "REUSED_CPIC" &
  !paste(NatMed_DrugCategory, NatMed_CategoryCode, sep = "||") %in%
    paste(
      known$NatMed_DrugCategory,
      known$NatMed_CategoryCode,
      sep = "||"
    )
]

if (nrow(bad)) {
  stop(
    "A reused CPIC mapping points to a category/code not present ",
    "in Gene_DrugModules.tsv.",
    call. = FALSE
  )
}

# Clean helper columns.
mapping[, c(
  "CPIC_MappingStatus",
  "CPIC_MappingBasis",
  "CPIC_MappingNotes"
) := NULL]

setcolorder(
  mapping,
  c(
    "DGIdb_Drug",
    "DGIdb_Drug_Normalized",
    "MappingStatus",
    "ReviewPriority",
    "NatMed_DrugCategory",
    "NatMed_CategoryCode",
    "CandidateNatMedHint",
    "MappingBasis",
    "MappingNotes",
    "CPIC_SourceDrug",
    "DGIdb_InteractionRows",
    "DGIdb_UniqueGenes",
    "DGIdb_InteractionTypes",
    "DGIdb_Sources",
    "DGIdb_MaxEvidenceScore",
    "DGIdb_MaxInteractionScore"
  )
)

priority_order <- c("ACCEPTED" = 1L, "REVIEW" = 2L, "LOW" = 3L)
mapping[, .priority_order := priority_order[ReviewPriority]]
setorder(
  mapping,
  .priority_order,
  NatMed_DrugCategory,
  CandidateNatMedHint,
  DGIdb_Drug
)
mapping[, .priority_order := NULL]

fwrite(
  mapping,
  out_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

qc <- data.table(
  Metric = c(
    "Raw DGIdb interaction rows",
    "Raw DGIdb unique drugs",
    "Raw DGIdb unique genes",
    "Approved interaction rows",
    "Approved unique drugs",
    "Approved unique genes",
    "Approved drugs reusing accepted CPIC mapping",
    "Additional DGIdb drugs flagged for NatMed review",
    "Unreviewed low-priority approved DGIdb drugs",
    "NatMed categories represented by reused mappings"
  ),
  Value = c(
    nrow(dg),
    uniqueN(dg$DrugNormalized),
    uniqueN(dg$gene_name),
    nrow(approved),
    uniqueN(approved$DrugNormalized),
    uniqueN(approved$gene_name),
    mapping[ReviewPriority == "ACCEPTED", .N],
    mapping[ReviewPriority == "REVIEW", .N],
    mapping[ReviewPriority == "LOW", .N],
    uniqueN(
      mapping[
        ReviewPriority == "ACCEPTED",
        NatMed_DrugCategory
      ]
    )
  )
)

fwrite(qc, qc_file, sep = "\t", quote = FALSE)

cat("\n============================================================\n")
cat("DGIdb mapping preparation completed successfully\n")
cat("============================================================\n\n")
print(qc)

cat("\nAdditional drugs requiring manual review:\n")
print(
  mapping[
    ReviewPriority == "REVIEW",
    .(
      DGIdb_Drug,
      CandidateNatMedHint,
      DGIdb_UniqueGenes,
      DGIdb_InteractionTypes,
      DGIdb_Sources
    )
  ]
)

cat(
  "\nOutputs:\n",
  "  ", out_file, "\n",
  "  ", qc_file, "\n\n",
  sep = ""
)

cat(
  "NEXT STEP:\n",
  "Review only rows with ReviewPriority == REVIEW.\n",
  "Do not run the DGIdb gene-drug benchmark until these mappings\n",
  "have been manually accepted/rejected.\n"
)
