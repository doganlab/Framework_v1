#!/usr/bin/env Rscript

# ============================================================
# Framework_v1
# 11_PGx_benchmark.R
#
# External validation of gene-drug modules against CPIC.
#
# Primary benchmark:
#   CPIC Level A/B gene-drug pairs whose individual drug maps
#   conservatively (EXACT or CLASS) to one of the observed
#   Nature Medicine drug categories.
#
# Secondary benchmark:
#   All accepted (EXACT/CLASS) CPIC pairs, with guideline,
#   used-for-recommendation and provisional fields retained.
#
# IMPORTANT
# ---------
# Recovery is evaluated at the NatMed drug-CATEGORY level.
# Example:
#   CPIC: SLCO1B1 - atorvastatin
#   NatMed: SLCO1B1 - lipophilic statin
# This is category-level recovery, not exact individual-drug
# recovery.
#
# No composite prioritisation score is calculated.
# ============================================================

options(stringsAsFactors = FALSE, warn = 1, scipen = 999)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
source(file.path(dirname(script_path), "..", "R", "bootstrap.R"))

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("The data.table package is required.", call. = FALSE)
}
library(data.table)

cpic_file <- file.path(DIR_PROCESSED, "CPIC_PGx_Benchmark.tsv")
mapping_file <- file.path(
  DIR_PROCESSED,
  "CPIC_to_NatMed_DrugCategory_Mapping.tsv"
)
module_file <- file.path(
  PROJECT_ROOT,
  "results",
  "gene_drug_modules",
  "Gene_DrugModules.tsv"
)

output_dir <- file.path(PROJECT_ROOT, "results", "pgx_benchmark")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

master_file <- file.path(output_dir, "PGxBenchmark_Master.tsv")
qc_file <- file.path(output_dir, "PGxBenchmark_QC.tsv")
collapsed_file <- file.path(output_dir, "PGxBenchmark_CategoryCollapsed.tsv")

required_files <- c(cpic_file, mapping_file, module_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop(
    "Missing required file(s):\n",
    paste(missing_files, collapse = "\n"),
    call. = FALSE
  )
}

cat("\n=== CPIC external PGx benchmark ===\n")

cpic <- fread(cpic_file, na.strings = c("", "NA"))
mapping <- fread(mapping_file, na.strings = c("", "NA"))
modules <- fread(module_file, na.strings = c("", "NA"))

required_cpic <- c(
  "Gene", "CPIC_Drug", "CPIC_Drug_Normalized",
  "CPIC_Level", "CPIC_Guideline",
  "HasCPICGuideline", "CPIC_A_or_B"
)

required_mapping <- c(
  "CPIC_Drug", "CPIC_Drug_Normalized",
  "NatMed_DrugCategory", "NatMed_CategoryCode",
  "MappingStatus", "MappingBasis", "MappingNotes"
)

required_modules <- c(
  "AssignedGene", "DrugCategory", "CategoryCode",
  "UniqueSNPs", "UniqueMetabolites",
  "Tier1SNPs", "Tier1Metabolites",
  "Tier2SNPs", "Tier2Metabolites",
  "CodingSNPs", "SpliceRelatedSNPs",
  "HighImpactSNPs", "ModerateImpactSNPs",
  "MinMRFDR", "MaxPPH4", "MinM2P",
  "SNPs", "Metabolites",
  "CardiovascularATCGroup"
)

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

check_cols(cpic, required_cpic, "CPIC_PGx_Benchmark.tsv")
check_cols(mapping, required_mapping,
           "CPIC_to_NatMed_DrugCategory_Mapping.tsv")
check_cols(modules, required_modules, "Gene_DrugModules.tsv")

# Normalize identifiers defensively.
cpic[, Gene := toupper(trimws(as.character(Gene)))]
cpic[, CPIC_Drug_Normalized :=
       tolower(trimws(as.character(CPIC_Drug_Normalized)))]

mapping[, CPIC_Drug_Normalized :=
          tolower(trimws(as.character(CPIC_Drug_Normalized)))]

modules[, AssignedGene :=
          toupper(trimws(as.character(AssignedGene)))]

# The reviewed mapping should have exactly one row per CPIC drug.
if (anyDuplicated(mapping$CPIC_Drug_Normalized) > 0) {
  dup <- unique(
    mapping[
      duplicated(CPIC_Drug_Normalized) |
        duplicated(CPIC_Drug_Normalized, fromLast = TRUE),
      CPIC_Drug_Normalized
    ]
  )
  stop(
    "Mapping file contains duplicated CPIC drugs: ",
    paste(head(dup, 20), collapse = ", "),
    call. = FALSE
  )
}

# Only manually accepted mappings are testable.
accepted_mapping <- mapping[
  MappingStatus %in% c("EXACT", "CLASS") &
    !is.na(NatMed_DrugCategory) &
    NatMed_DrugCategory != ""
]

if (nrow(accepted_mapping) == 0) {
  stop(
    "No EXACT/CLASS mappings found. Did the reviewed mapping replace ",
    "data/processed/CPIC_to_NatMed_DrugCategory_Mapping.tsv?",
    call. = FALSE
  )
}

# Add mapping to CPIC pairs.
bench <- merge(
  cpic,
  accepted_mapping[, .(
    CPIC_Drug_Normalized,
    NatMed_DrugCategory,
    NatMed_CategoryCode,
    MappingStatus,
    MappingBasis,
    MappingNotes
  )],
  by = "CPIC_Drug_Normalized",
  all = FALSE,
  sort = FALSE
)

# One CPIC gene-drug pair should remain one row.
if (anyDuplicated(
  bench[, paste(Gene, CPIC_Drug_Normalized, sep = "||")]
) > 0) {
  stop(
    "Duplicate CPIC gene-drug pairs appeared after mapping.",
    call. = FALSE
  )
}

# Confirm accepted mapping points to categories that actually exist
# in the gene-drug module table.
known_categories <- unique(
  modules[, .(DrugCategory, CategoryCode)]
)

bad_category <- bench[
  !paste(NatMed_DrugCategory, NatMed_CategoryCode, sep = "||") %in%
    paste(known_categories$DrugCategory,
          known_categories$CategoryCode, sep = "||")
]

if (nrow(bad_category) > 0) {
  stop(
    "Accepted mapping contains NatMed category/code combinations ",
    "not present in Gene_DrugModules.tsv. First examples:\n",
    paste(
      head(unique(
        paste(
          bad_category$NatMed_DrugCategory,
          bad_category$NatMed_CategoryCode,
          sep = " | "
        )
      ), 10),
      collapse = "\n"
    ),
    call. = FALSE
  )
}

# Prepare one module row per gene/category.
module_key <- modules[, .(
  Gene = AssignedGene,
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
  module_key[, paste(Gene, NatMed_DrugCategory,
                     NatMed_CategoryCode, sep = "||")]
) > 0) {
  stop(
    "Gene_DrugModules.tsv contains duplicate gene-category modules.",
    call. = FALSE
  )
}

# Left join: every testable CPIC pair remains in the denominator.
master <- merge(
  bench,
  module_key,
  by = c("Gene", "NatMed_DrugCategory", "NatMed_CategoryCode"),
  all.x = TRUE,
  sort = FALSE
)

master[, Recovered := !is.na(UniqueSNPs)]

# Explicit benchmark strata.
master[, Primary_AB := CPIC_Level %in% c("A", "B")]
master[, GuidelineLinked := HasCPICGuideline %in% TRUE]

# Retain extra CPIC fields when present.
preferred_front <- c(
  "Gene",
  "CPIC_Drug",
  "CPIC_Drug_Normalized",
  "CPIC_Level",
  "CPIC_Guideline",
  "HasCPICGuideline",
  "CPIC_A_or_B",
  "Primary_AB",
  "NatMed_DrugCategory",
  "NatMed_CategoryCode",
  "MappingStatus",
  "MappingBasis",
  "MappingNotes",
  "Recovered",
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
)

preferred_front <- intersect(preferred_front, names(master))
remaining <- setdiff(names(master), preferred_front)
setcolorder(master, c(preferred_front, remaining))

setorder(
  master,
  -Primary_AB,
  -Recovered,
  CPIC_Level,
  Gene,
  CPIC_Drug
)

fwrite(
  master,
  master_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# ------------------------------------------------------------
# NatMed-resolution category-collapsed benchmark
# ------------------------------------------------------------
# Multiple CPIC individual drugs can map to the same NatMed category
# (e.g. several statins). Collapse them to one Gene + NatMed category
# unit so that the primary validation is evaluated at the same
# resolution as the pharmacometabolomic evidence.

collapsed <- master[, .(
  CPIC_Drugs = paste(sort(unique(CPIC_Drug)), collapse = ";"),
  N_CPIC_Drugs = uniqueN(CPIC_Drug_Normalized),
  CPIC_Levels = paste(sort(unique(na.omit(CPIC_Level))), collapse = ";"),
  Any_Level_A = isTRUE(any(CPIC_Level == "A", na.rm = TRUE)),
  Any_Level_B = isTRUE(any(CPIC_Level == "B", na.rm = TRUE)),
  Primary_AB = isTRUE(any(CPIC_Level %in% c("A", "B"), na.rm = TRUE)),
  GuidelineLinked = isTRUE(any(GuidelineLinked == TRUE, na.rm = TRUE)),
  MappingStatuses = paste(sort(unique(MappingStatus)), collapse = ";"),
  Recovered = isTRUE(any(Recovered == TRUE, na.rm = TRUE)),

  # Module evidence is identical within a recovered Gene-category
  # combination; use first non-missing value defensively.
  UniqueSNPs = {
    z <- na.omit(UniqueSNPs); if (length(z)) as.numeric(z[1]) else NA_real_
  },
  UniqueMetabolites = {
    z <- na.omit(UniqueMetabolites); if (length(z)) as.numeric(z[1]) else NA_real_
  },
  Tier1SNPs = {
    z <- na.omit(Tier1SNPs); if (length(z)) as.numeric(z[1]) else NA_real_
  },
  Tier1Metabolites = {
    z <- na.omit(Tier1Metabolites); if (length(z)) as.numeric(z[1]) else NA_real_
  },
  CodingSNPs = {
    z <- na.omit(CodingSNPs); if (length(z)) as.numeric(z[1]) else NA_real_
  },
  SpliceRelatedSNPs = {
    z <- na.omit(SpliceRelatedSNPs); if (length(z)) as.numeric(z[1]) else NA_real_
  },
  HighImpactSNPs = {
    z <- na.omit(HighImpactSNPs); if (length(z)) as.numeric(z[1]) else NA_real_
  },
  MinMRFDR = {
    z <- na.omit(MinMRFDR); if (length(z)) as.numeric(z[1]) else NA_real_
  },
  MaxPPH4 = {
    z <- na.omit(MaxPPH4); if (length(z)) as.numeric(z[1]) else NA_real_
  },
  MinM2P = {
    z <- na.omit(MinM2P); if (length(z)) as.numeric(z[1]) else NA_real_
  },
  SNPs = {
    z <- na.omit(SNPs); if (length(z)) z[1] else NA_character_
  },
  Metabolites = {
    z <- na.omit(Metabolites); if (length(z)) z[1] else NA_character_
  }
), by = .(Gene, NatMed_DrugCategory, NatMed_CategoryCode)]

setorder(collapsed, -Primary_AB, -Recovered, Gene, NatMed_DrugCategory)

fwrite(
  collapsed,
  collapsed_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# ------------------------------------------------------------
# QC / recovery metrics
# ------------------------------------------------------------

safe_pct <- function(n, d) {
  if (d == 0) return(NA_real_)
  100 * n / d
}

all_n <- nrow(master)
all_rec <- master[Recovered == TRUE, .N]

ab_n <- master[Primary_AB == TRUE, .N]
ab_rec <- master[Primary_AB == TRUE & Recovered == TRUE, .N]

guideline_n <- master[GuidelineLinked == TRUE, .N]
guideline_rec <- master[
  GuidelineLinked == TRUE & Recovered == TRUE, .N
]

a_n <- master[CPIC_Level == "A", .N]
a_rec <- master[CPIC_Level == "A" & Recovered == TRUE, .N]

b_n <- master[CPIC_Level == "B", .N]
b_rec <- master[CPIC_Level == "B" & Recovered == TRUE, .N]

collapsed_n <- nrow(collapsed)
collapsed_rec <- collapsed[Recovered == TRUE, .N]
collapsed_ab_n <- collapsed[Primary_AB == TRUE, .N]
collapsed_ab_rec <- collapsed[
  Primary_AB == TRUE & Recovered == TRUE, .N
]
collapsed_a_n <- collapsed[Any_Level_A == TRUE, .N]
collapsed_a_rec <- collapsed[
  Any_Level_A == TRUE & Recovered == TRUE, .N
]

exact_n <- master[MappingStatus == "EXACT", .N]
exact_rec <- master[
  MappingStatus == "EXACT" & Recovered == TRUE, .N
]

class_n <- master[MappingStatus == "CLASS", .N]
class_rec <- master[
  MappingStatus == "CLASS" & Recovered == TRUE, .N
]

qc <- data.table(
  Metric = c(
    "Accepted CPIC drugs mapped to NatMed",
    "NatMed categories represented in CPIC benchmark",
    "All testable CPIC gene-drug pairs",
    "All testable CPIC pairs recovered",
    "All-pair recovery percent",
    "Primary CPIC Level A/B testable pairs",
    "Primary CPIC Level A/B pairs recovered",
    "Primary A/B recovery percent",
    "CPIC Level A testable pairs",
    "CPIC Level A pairs recovered",
    "CPIC Level A recovery percent",
    "CPIC Level B testable pairs",
    "CPIC Level B pairs recovered",
    "CPIC Level B recovery percent",
    "Category-collapsed testable gene-drug modules",
    "Category-collapsed modules recovered",
    "Category-collapsed recovery percent",
    "Category-collapsed Level A/B modules",
    "Category-collapsed Level A/B modules recovered",
    "Category-collapsed Level A/B recovery percent",
    "Category-collapsed Level A modules",
    "Category-collapsed Level A modules recovered",
    "Category-collapsed Level A recovery percent",
    "Guideline-linked testable pairs",
    "Guideline-linked pairs recovered",
    "Guideline-linked recovery percent",
    "EXACT-mapped testable pairs",
    "EXACT-mapped pairs recovered",
    "CLASS-mapped testable pairs",
    "CLASS-mapped pairs recovered",
    "Recovered pairs with coding SNP support",
    "Recovered pairs with splice-related SNP support",
    "Recovered pairs with HIGH-impact SNP support"
  ),
  Value = c(
    as.character(uniqueN(master$CPIC_Drug_Normalized)),
    as.character(uniqueN(master$NatMed_DrugCategory)),
    as.character(all_n),
    as.character(all_rec),
    sprintf("%.2f", safe_pct(all_rec, all_n)),
    as.character(ab_n),
    as.character(ab_rec),
    sprintf("%.2f", safe_pct(ab_rec, ab_n)),
    as.character(a_n),
    as.character(a_rec),
    sprintf("%.2f", safe_pct(a_rec, a_n)),
    as.character(b_n),
    as.character(b_rec),
    sprintf("%.2f", safe_pct(b_rec, b_n)),
    as.character(collapsed_n),
    as.character(collapsed_rec),
    sprintf("%.2f", safe_pct(collapsed_rec, collapsed_n)),
    as.character(collapsed_ab_n),
    as.character(collapsed_ab_rec),
    sprintf("%.2f", safe_pct(collapsed_ab_rec, collapsed_ab_n)),
    as.character(collapsed_a_n),
    as.character(collapsed_a_rec),
    sprintf("%.2f", safe_pct(collapsed_a_rec, collapsed_a_n)),
    as.character(guideline_n),
    as.character(guideline_rec),
    sprintf("%.2f", safe_pct(guideline_rec, guideline_n)),
    as.character(exact_n),
    as.character(exact_rec),
    as.character(class_n),
    as.character(class_rec),
    as.character(
      master[Recovered == TRUE & CodingSNPs > 0, .N]
    ),
    as.character(
      master[Recovered == TRUE & SpliceRelatedSNPs > 0, .N]
    ),
    as.character(
      master[Recovered == TRUE & HighImpactSNPs > 0, .N]
    )
  )
)

fwrite(qc, qc_file, sep = "\t", quote = FALSE)

cat("\n============================================================\n")
cat("CPIC PGx benchmark completed successfully\n")
cat("============================================================\n\n")
print(qc)

cat("\nPrimary Level A/B benchmark:\n")
cat(
  "  Recovered ",
  ab_rec, " / ", ab_n,
  " (", sprintf("%.1f", safe_pct(ab_rec, ab_n)), "%)\n",
  sep = ""
)

cat("\nCategory-collapsed primary Level A/B benchmark:\n")
cat(
  "  Recovered ",
  collapsed_ab_rec, " / ", collapsed_ab_n,
  " (", sprintf("%.1f", safe_pct(collapsed_ab_rec, collapsed_ab_n)), "%)\n",
  sep = ""
)

cat("\nCategory-collapsed Level A/B module details:\n")
print(
  collapsed[
    Primary_AB == TRUE,
    .(
      Gene,
      CPIC_Drugs,
      CPIC_Levels,
      NatMed_DrugCategory,
      Recovered,
      UniqueSNPs,
      UniqueMetabolites,
      CodingSNPs,
      SpliceRelatedSNPs
    )
  ]
)

cat("\nIndividual Level A/B pair details (supplementary resolution):\n")
print(
  master[
    Primary_AB == TRUE,
    .(
      Gene,
      CPIC_Drug,
      CPIC_Level,
      NatMed_DrugCategory,
      MappingStatus,
      Recovered,
      UniqueSNPs,
      UniqueMetabolites,
      CodingSNPs,
      SpliceRelatedSNPs
    )
  ]
)

cat(
  "\nOutputs:\n",
  "  ", master_file, "\n",
  "  ", collapsed_file, "\n",
  "  ", qc_file, "\n\n",
  sep = ""
)

cat(
  "INTERPRETATION NOTE:\n",
  "Recovery is at the NatMed drug-category level. A recovered CPIC\n",
  "individual-drug pair therefore indicates that the same gene was\n",
  "recovered for the corresponding NatMed drug category; it is not\n",
  "evidence that the framework resolved the individual drug itself.\n"
)
