#!/usr/bin/env Rscript

# ============================================================
# Framework_v1
# 15_ChEMBL_benchmark.R
#
# External direct mechanism-of-action validation against ChEMBL.
#
# PRIMARY benchmark:
#   ChEMBL SINGLE PROTEIN target + NatMed drug category
#
# EXTENDED benchmark:
#   ChEMBL target + NatMed drug category
#   (SINGLE PROTEIN, PROTEIN COMPLEX, PROTEIN FAMILY, etc.)
#
# For complex/family targets, recovery is assessed at the TARGET
# level: a target is recovered if at least one of its annotated
# human component genes is present in our Gene_DrugModules for the
# same NatMed drug category. A 40-component complex therefore counts
# as ONE ChEMBL target-category validation unit, not 40 successes.
#
# No composite score is calculated.
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

chembl_evidence_file <- file.path(
  DIR_PROCESSED,
  "ChEMBL_MoA_GeneDrugEvidence.tsv"
)

module_file <- file.path(
  PROJECT_ROOT,
  "results",
  "gene_drug_modules",
  "Gene_DrugModules.tsv"
)

output_dir <- file.path(
  PROJECT_ROOT,
  "results",
  "chembl_benchmark"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

target_master_file <- file.path(
  output_dir,
  "ChEMBLBenchmark_TargetMaster.tsv"
)

gene_supp_file <- file.path(
  output_dir,
  "ChEMBLBenchmark_GeneSupplement.tsv"
)

qc_file <- file.path(
  output_dir,
  "ChEMBLBenchmark_QC.tsv"
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

safe_pct <- function(n, d) {
  if (d == 0) return(NA_real_)
  100 * n / d
}

first_nonmissing_num <- function(x) {
  z <- suppressWarnings(as.numeric(x))
  z <- z[!is.na(z)]
  if (!length(z)) NA_real_ else z[1]
}

# ------------------------------------------------------------
# Validate and load
# ------------------------------------------------------------

required_files <- c(chembl_evidence_file, module_file)
missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing required file(s):\n",
    paste(missing_files, collapse = "\n"),
    call. = FALSE
  )
}

cat("\n=== ChEMBL direct mechanism-of-action benchmark ===\n")

chembl <- fread(
  chembl_evidence_file,
  sep = "\t",
  na.strings = c("", "NA"),
  showProgress = FALSE
)

modules <- fread(
  module_file,
  sep = "\t",
  na.strings = c("", "NA"),
  showProgress = FALSE
)

check_cols(
  chembl,
  c(
    "Gene",
    "InputDrug",
    "NatMed_DrugCategory",
    "NatMed_CategoryCode",
    "MoleculeChEMBLID",
    "TargetChEMBLID",
    "TargetPrefName",
    "TargetType",
    "ActionType",
    "MechanismOfAction",
    "DirectInteraction",
    "MolecularMechanism",
    "DiseaseEfficacy",
    "Accession",
    "MechanismRefs"
  ),
  "ChEMBL_MoA_GeneDrugEvidence.tsv"
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
    "Metabolites"
  ),
  "Gene_DrugModules.tsv"
)

chembl[, Gene := toupper(trimws(as.character(Gene)))]
modules[, AssignedGene :=
          toupper(trimws(as.character(AssignedGene)))]

chembl <- chembl[
  !is.na(Gene) & Gene != "" &
  !is.na(TargetChEMBLID) & TargetChEMBLID != "" &
  !is.na(NatMed_DrugCategory) & NatMed_DrugCategory != ""
]

# ------------------------------------------------------------
# Prepare our module lookup
# ------------------------------------------------------------

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
  Metabolites
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
    "Gene_DrugModules.tsv contains duplicate gene-category modules.",
    call. = FALSE
  )
}

# ------------------------------------------------------------
# 1. Supplementary component-gene benchmark
# ------------------------------------------------------------
# This table is useful for mechanistic inspection, but is NOT the
# primary denominator for complex/family targets.

gene_units <- unique(
  chembl[, .(
    Gene,
    NatMed_DrugCategory,
    NatMed_CategoryCode,
    TargetChEMBLID,
    TargetPrefName,
    TargetType,
    ActionType,
    MechanismOfAction,
    InputDrug,
    MoleculeChEMBLID,
    Accession,
    DirectInteraction,
    MolecularMechanism,
    DiseaseEfficacy,
    MechanismRefs
  )]
)

gene_supp <- merge(
  gene_units,
  module_key,
  by = c(
    "Gene",
    "NatMed_DrugCategory",
    "NatMed_CategoryCode"
  ),
  all.x = TRUE,
  sort = FALSE
)

gene_supp[, Recovered := !is.na(UniqueSNPs)]
gene_supp[, PrimarySingleProtein :=
            TargetType == "SINGLE PROTEIN"]

setorder(
  gene_supp,
  -PrimarySingleProtein,
  -Recovered,
  Gene,
  NatMed_DrugCategory,
  TargetChEMBLID
)

fwrite(
  gene_supp,
  gene_supp_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# ------------------------------------------------------------
# 2. Target-level benchmark
# ------------------------------------------------------------
# Collapse all component genes belonging to the same ChEMBL target
# and NatMed category. This prevents large protein complexes from
# inflating the validation denominator.

target_units <- chembl[
  ,
  .(
    ChEMBL_ComponentGenes = collapse_unique(Gene, 100),
    N_ComponentGenes = uniqueN(Gene),
    ChEMBL_InputDrugs = collapse_unique(InputDrug, 50),
    N_InputDrugs = uniqueN(InputDrug),
    ChEMBL_Molecules = collapse_unique(MoleculeChEMBLID, 50),
    N_ChEMBLMolecules = uniqueN(MoleculeChEMBLID),
    ActionTypes = collapse_unique(ActionType, 20),
    MechanismsOfAction = collapse_unique(MechanismOfAction, 20),
    Accessions = collapse_unique(Accession, 100),
    AnyDirectInteraction = any(
      as.character(DirectInteraction) %in% c("1", "TRUE", "true"),
      na.rm = TRUE
    ),
    AnyMolecularMechanism = any(
      as.character(MolecularMechanism) %in% c("1", "TRUE", "true"),
      na.rm = TRUE
    ),
    AnyDiseaseEfficacy = any(
      as.character(DiseaseEfficacy) %in% c("1", "TRUE", "true"),
      na.rm = TRUE
    ),
    MechanismRefs = collapse_unique(MechanismRefs, 20)
  ),
  by = .(
    TargetChEMBLID,
    TargetPrefName,
    TargetType,
    NatMed_DrugCategory,
    NatMed_CategoryCode
  )
]

if (anyDuplicated(
  target_units[
    ,
    paste(
      TargetChEMBLID,
      NatMed_DrugCategory,
      NatMed_CategoryCode,
      sep = "||"
    )
  ]
) > 0) {
  stop(
    "Duplicate ChEMBL target-category units remain after collapse.",
    call. = FALSE
  )
}

# ------------------------------------------------------------
# Determine recovered component genes per target-category unit.
# ------------------------------------------------------------

target_recovery <- gene_supp[
  ,
  .(
    Recovered = any(Recovered == TRUE, na.rm = TRUE),
    RecoveredGenes = collapse_unique(
      Gene[Recovered == TRUE],
      100
    ),
    N_RecoveredGenes = uniqueN(
      Gene[Recovered == TRUE]
    ),

    # Aggregate our evidence only over RECOVERED component genes.
    TotalRecoveredGeneSNPs = sum(
      UniqueSNPs[Recovered == TRUE],
      na.rm = TRUE
    ),
    MaxUniqueSNPsPerRecoveredGene = {
      z <- suppressWarnings(
        as.numeric(UniqueSNPs[Recovered == TRUE])
      )
      z <- z[!is.na(z)]
      if (!length(z)) NA_real_ else max(z)
    },
    MaxUniqueMetabolitesPerRecoveredGene = {
      z <- suppressWarnings(
        as.numeric(UniqueMetabolites[Recovered == TRUE])
      )
      z <- z[!is.na(z)]
      if (!length(z)) NA_real_ else max(z)
    },
    RecoveredGenesWithCodingSNP = uniqueN(
      Gene[
        Recovered == TRUE &
          !is.na(CodingSNPs) &
          CodingSNPs > 0
      ]
    ),
    RecoveredGenesWithSpliceSNP = uniqueN(
      Gene[
        Recovered == TRUE &
          !is.na(SpliceRelatedSNPs) &
          SpliceRelatedSNPs > 0
      ]
    ),
    RecoveredGenesWithHighImpactSNP = uniqueN(
      Gene[
        Recovered == TRUE &
          !is.na(HighImpactSNPs) &
          HighImpactSNPs > 0
      ]
    ),
    RecoveredSNPs = collapse_unique(
      SNPs[Recovered == TRUE],
      50
    ),
    RecoveredMetabolites = collapse_unique(
      Metabolites[Recovered == TRUE],
      50
    )
  ),
  by = .(
    TargetChEMBLID,
    NatMed_DrugCategory,
    NatMed_CategoryCode
  )
]

target_master <- merge(
  target_units,
  target_recovery,
  by = c(
    "TargetChEMBLID",
    "NatMed_DrugCategory",
    "NatMed_CategoryCode"
  ),
  all.x = TRUE,
  sort = FALSE
)

target_master[is.na(Recovered), Recovered := FALSE]
target_master[is.na(N_RecoveredGenes), N_RecoveredGenes := 0L]

target_master[, BenchmarkStratum := fifelse(
  TargetType == "SINGLE PROTEIN",
  "Primary_SingleProtein",
  "Extended_ComplexFamilyOther"
)]

setorder(
  target_master,
  BenchmarkStratum,
  -Recovered,
  -N_RecoveredGenes,
  NatMed_DrugCategory,
  TargetPrefName
)

fwrite(
  target_master,
  target_master_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# ------------------------------------------------------------
# QC / benchmark summary
# ------------------------------------------------------------

primary <- target_master[
  BenchmarkStratum == "Primary_SingleProtein"
]

extended_non_single <- target_master[
  BenchmarkStratum == "Extended_ComplexFamilyOther"
]

all_targets <- target_master

primary_n <- nrow(primary)
primary_rec <- primary[Recovered == TRUE, .N]

extended_n <- nrow(all_targets)
extended_rec <- all_targets[Recovered == TRUE, .N]

non_single_n <- nrow(extended_non_single)
non_single_rec <- extended_non_single[
  Recovered == TRUE,
  .N
]

qc <- data.table(
  Metric = c(
    "ChEMBL component-gene evidence rows",
    "Unique human component genes",
    "Unique ChEMBL mechanism targets",
    "Unique NatMed categories in ChEMBL MoA",
    "Primary SINGLE PROTEIN target-category units",
    "Primary SINGLE PROTEIN units recovered",
    "Primary SINGLE PROTEIN recovery percent",
    "Extended all-target category units",
    "Extended all-target units recovered",
    "Extended all-target recovery percent",
    "Non-single-protein target-category units",
    "Non-single-protein units recovered",
    "Recovered target units with coding-SNP gene support",
    "Recovered target units with splice-SNP gene support",
    "Recovered target units with HIGH-impact-SNP gene support"
  ),
  Value = c(
    as.character(nrow(gene_supp)),
    as.character(uniqueN(gene_supp$Gene)),
    as.character(uniqueN(target_master$TargetChEMBLID)),
    as.character(uniqueN(target_master$NatMed_DrugCategory)),
    as.character(primary_n),
    as.character(primary_rec),
    sprintf("%.2f", safe_pct(primary_rec, primary_n)),
    as.character(extended_n),
    as.character(extended_rec),
    sprintf("%.2f", safe_pct(extended_rec, extended_n)),
    as.character(non_single_n),
    as.character(non_single_rec),
    as.character(
      target_master[
        Recovered == TRUE &
          RecoveredGenesWithCodingSNP > 0,
        .N
      ]
    ),
    as.character(
      target_master[
        Recovered == TRUE &
          RecoveredGenesWithSpliceSNP > 0,
        .N
      ]
    ),
    as.character(
      target_master[
        Recovered == TRUE &
          RecoveredGenesWithHighImpactSNP > 0,
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
cat("ChEMBL MoA benchmark completed successfully\n")
cat("============================================================\n\n")
print(qc)

cat("\nPrimary SINGLE PROTEIN benchmark:\n")
cat(
  "  Recovered ",
  primary_rec, " / ", primary_n,
  " (",
  sprintf("%.1f", safe_pct(primary_rec, primary_n)),
  "%)\n",
  sep = ""
)

cat("\nPrimary SINGLE PROTEIN target details:\n")
print(
  primary[
    ,
    .(
      TargetPrefName,
      ChEMBL_ComponentGenes,
      NatMed_DrugCategory,
      ActionTypes,
      MechanismsOfAction,
      Recovered,
      RecoveredGenes,
      N_RecoveredGenes,
      MaxUniqueSNPsPerRecoveredGene,
      MaxUniqueMetabolitesPerRecoveredGene
    )
  ]
)

cat("\nExtended non-single-protein target details:\n")
print(
  extended_non_single[
    ,
    .(
      TargetPrefName,
      TargetType,
      N_ComponentGenes,
      NatMed_DrugCategory,
      ActionTypes,
      Recovered,
      RecoveredGenes,
      N_RecoveredGenes
    )
  ]
)

cat(
  "\nOutputs:\n",
  "  ", target_master_file, "\n",
  "  ", gene_supp_file, "\n",
  "  ", qc_file, "\n\n",
  sep = ""
)

cat(
  "INTERPRETATION NOTE:\n",
  "Primary validation uses SINGLE PROTEIN targets only. Complex/family\n",
  "targets are counted once per ChEMBL target + NatMed category in the\n",
  "extended analysis, regardless of the number of component genes.\n",
  "ChEMBL MoA recovery is direct pharmacological mechanism support and\n",
  "must not be interpreted as clinical PGx actionability.\n"
)
