#!/usr/bin/env Rscript

# ============================================================
# Framework_v1
# 16_integrated_evidence_summary.R
#
# Publication-oriented integration of the final gene-drug modules
# with the three independent external validation layers:
#
#   CPIC   -> clinical pharmacogenetic support
#   DGIdb  -> broad pharmacological gene-drug interaction support
#   ChEMBL -> curated direct mechanism-of-action support
#
# The base analytical unit is:
#
#   Gene + NatMed drug category
#
# IMPORTANT
# ---------
# No composite prioritisation score is calculated.
#
# ExternalDatabaseOverlapCount is a DESCRIPTIVE count of database overlaps
# (CPIC-tracked, DGIdb, ChEMBL), not a ranking score.
#
# The internal genetic/metabolomic evidence dimensions remain in
# separate columns (SNPs, metabolites, Tier1 support, coding/splice,
# MR/colocalization/drug-metabolite evidence).
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

module_file <- file.path(
  PROJECT_ROOT,
  "results",
  "gene_drug_modules",
  "Gene_DrugModules.tsv"
)

cpic_file <- file.path(
  PROJECT_ROOT,
  "results",
  "pgx_benchmark",
  "PGxBenchmark_CategoryCollapsed.tsv"
)

dgidb_file <- file.path(
  PROJECT_ROOT,
  "results",
  "dgidb_benchmark",
  "DGIdbBenchmark_Master.tsv"
)

chembl_gene_file <- file.path(
  PROJECT_ROOT,
  "results",
  "chembl_benchmark",
  "ChEMBLBenchmark_GeneSupplement.tsv"
)

output_dir <- file.path(
  PROJECT_ROOT,
  "results",
  "integrated_evidence"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

master_file <- file.path(
  output_dir,
  "IntegratedGeneDrugEvidence.tsv"
)

qc_file <- file.path(
  output_dir,
  "IntegratedGeneDrugEvidence_QC.tsv"
)

pattern_file <- file.path(
  output_dir,
  "ExternalDatabasePattern_Summary.tsv"
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

flag01 <- function(x) as.integer(x %in% TRUE)

# ------------------------------------------------------------
# Validate input files
# ------------------------------------------------------------

required_files <- c(
  module_file,
  cpic_file,
  dgidb_file,
  chembl_gene_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing required file(s):\n",
    paste(missing_files, collapse = "\n"),
    call. = FALSE
  )
}

cat("\n=== Integrated external evidence summary ===\n")

modules <- fread(
  module_file,
  sep = "\t",
  na.strings = c("", "NA"),
  showProgress = FALSE
)

cpic <- fread(
  cpic_file,
  sep = "\t",
  na.strings = c("", "NA"),
  showProgress = FALSE
)

dgidb <- fread(
  dgidb_file,
  sep = "\t",
  na.strings = c("", "NA"),
  showProgress = FALSE
)

chembl <- fread(
  chembl_gene_file,
  sep = "\t",
  na.strings = c("", "NA"),
  showProgress = FALSE
)

# ------------------------------------------------------------
# Validate required columns
# ------------------------------------------------------------

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
    "OverlappingGeneSNPs",
    "NearestGeneSNPs",
    "MinMRFDR",
    "MaxPPH4",
    "MinM2P",
    "SNPs",
    "Metabolites",
    "CardiovascularATCGroup"
  ),
  "Gene_DrugModules.tsv"
)

check_cols(
  cpic,
  c(
    "Gene",
    "NatMed_DrugCategory",
    "NatMed_CategoryCode",
    "CPIC_Drugs",
    "CPIC_Levels",
    "Primary_AB",
    "GuidelineLinked",
    "Recovered"
  ),
  "PGxBenchmark_CategoryCollapsed.tsv"
)

check_cols(
  dgidb,
  c(
    "Gene",
    "NatMed_DrugCategory",
    "NatMed_CategoryCode",
    "DGIdb_UniqueDrugs",
    "DGIdb_Drugs",
    "DGIdb_UniqueSources",
    "DGIdb_Sources",
    "DGIdb_InteractionTypes",
    "DGIdb_MaxEvidenceScore",
    "DGIdb_MaxInteractionScore",
    "Recovered"
  ),
  "DGIdbBenchmark_Master.tsv"
)

check_cols(
  chembl,
  c(
    "Gene",
    "NatMed_DrugCategory",
    "NatMed_CategoryCode",
    "TargetChEMBLID",
    "TargetPrefName",
    "TargetType",
    "ActionType",
    "MechanismOfAction",
    "InputDrug",
    "Recovered",
    "PrimarySingleProtein"
  ),
  "ChEMBLBenchmark_GeneSupplement.tsv"
)

# ------------------------------------------------------------
# Base table = our gene-drug modules
# ------------------------------------------------------------

base <- copy(modules)

setnames(
  base,
  old = c("AssignedGene", "DrugCategory", "CategoryCode"),
  new = c("Gene", "NatMed_DrugCategory", "NatMed_CategoryCode")
)

base[, Gene := toupper(trimws(as.character(Gene)))]

if (anyDuplicated(
  base[
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
    "Gene_DrugModules.tsv contains duplicate Gene + NatMed category rows.",
    call. = FALSE
  )
}

# ------------------------------------------------------------
# CPIC layer
# ------------------------------------------------------------
# Presence in the collapsed CPIC benchmark defines external CPIC
# support for the same Gene + NatMed category. Because BASE already
# contains only recovered framework modules, CPIC 'Recovered' itself
# is not used as the support definition after the join.

cpic_layer <- unique(
  cpic[
    ,
    .(
      Gene = toupper(trimws(as.character(Gene))),
      NatMed_DrugCategory,
      NatMed_CategoryCode,
      CPICTrackedPair = TRUE,
      CPIC_PrimaryAB = Primary_AB %in% TRUE,
      CPIC_GuidelineLinked = GuidelineLinked %in% TRUE,
      CPIC_Drugs,
      CPIC_Levels
    )
  ],
  by = c("Gene", "NatMed_DrugCategory", "NatMed_CategoryCode")
)

# ------------------------------------------------------------
# DGIdb layer
# ------------------------------------------------------------

dgidb_layer <- unique(
  dgidb[
    ,
    .(
      Gene = toupper(trimws(as.character(Gene))),
      NatMed_DrugCategory,
      NatMed_CategoryCode,
      DGIdb_Support = TRUE,
      DGIdb_UniqueDrugs,
      DGIdb_Drugs,
      DGIdb_UniqueSources,
      DGIdb_Sources,
      DGIdb_InteractionTypes,
      DGIdb_MaxEvidenceScore,
      DGIdb_MaxInteractionScore
    )
  ],
  by = c("Gene", "NatMed_DrugCategory", "NatMed_CategoryCode")
)

# ------------------------------------------------------------
# ChEMBL layer
# ------------------------------------------------------------
# Aggregate possible multiple ChEMBL target records per Gene+category.

chembl_layer <- chembl[
  ,
  .(
    ChEMBL_Support = TRUE,
    ChEMBL_PrimarySingleProteinSupport = any(
      PrimarySingleProtein %in% TRUE,
      na.rm = TRUE
    ),
    ChEMBL_Targets = collapse_unique(TargetChEMBLID, 20),
    ChEMBL_TargetNames = collapse_unique(TargetPrefName, 20),
    ChEMBL_TargetTypes = collapse_unique(TargetType, 20),
    ChEMBL_ActionTypes = collapse_unique(ActionType, 20),
    ChEMBL_Mechanisms = collapse_unique(MechanismOfAction, 20),
    ChEMBL_InputDrugs = collapse_unique(InputDrug, 30),
    ChEMBL_TargetCount = uniqueN(TargetChEMBLID)
  ),
  by = .(
    Gene = toupper(trimws(as.character(Gene))),
    NatMed_DrugCategory,
    NatMed_CategoryCode
  )
]

# ------------------------------------------------------------
# Integrate
# ------------------------------------------------------------

master <- merge(
  base,
  cpic_layer,
  by = c("Gene", "NatMed_DrugCategory", "NatMed_CategoryCode"),
  all.x = TRUE,
  sort = FALSE
)

master <- merge(
  master,
  dgidb_layer,
  by = c("Gene", "NatMed_DrugCategory", "NatMed_CategoryCode"),
  all.x = TRUE,
  sort = FALSE
)

master <- merge(
  master,
  chembl_layer,
  by = c("Gene", "NatMed_DrugCategory", "NatMed_CategoryCode"),
  all.x = TRUE,
  sort = FALSE
)

# Convert missing support flags to FALSE.
logical_flags <- c(
  "CPICTrackedPair",
  "CPIC_PrimaryAB",
  "CPIC_GuidelineLinked",
  "DGIdb_Support",
  "ChEMBL_Support",
  "ChEMBL_PrimarySingleProteinSupport"
)

for (nm in logical_flags) {
  if (nm %in% names(master)) {
    master[is.na(get(nm)), (nm) := FALSE]
  }
}

# ------------------------------------------------------------
# Descriptive evidence-axis fields
# ------------------------------------------------------------

master[, ExternalDatabaseOverlapCount :=
  flag01(CPICTrackedPair) +
  flag01(DGIdb_Support) +
  flag01(ChEMBL_Support)
]

master[, IndependentValidatedEvidenceCount :=
  flag01(CPIC_PrimaryAB) +
  flag01(DGIdb_Support) +
  flag01(ChEMBL_PrimarySingleProteinSupport)
]

master[, ExternalDatabasePattern := fifelse(
  ExternalDatabaseOverlapCount == 0,
  "None",
  paste0(
    fifelse(CPICTrackedPair, "CPIC", ""),
    fifelse(
      CPICTrackedPair & (DGIdb_Support | ChEMBL_Support),
      "+",
      ""
    ),
    fifelse(DGIdb_Support, "DGIdb", ""),
    fifelse(
      DGIdb_Support & ChEMBL_Support,
      "+",
      ""
    ),
    fifelse(ChEMBL_Support, "ChEMBL", "")
  )
)]

# More explicit interpretive flags.
master[, CPICTrackedPair := CPICTrackedPair]
master[, CPICLevelABSupport := CPIC_PrimaryAB]
master[, DGIdbSupport := DGIdb_Support]
master[, ChEMBLMoASupport := ChEMBL_Support]
master[, ChEMBLSingleProteinMoASupport :=
          ChEMBL_PrimarySingleProteinSupport]

master[, MultiMetaboliteInternalSupport :=
          UniqueMetabolites >= 2]

master[, MultiSNPInternalSupport :=
          UniqueSNPs >= 2]

master[, FunctionalVariantInternalSupport :=
          CodingSNPs > 0 | SpliceRelatedSNPs > 0]

# ------------------------------------------------------------
# Ordering only -- NOT a score-based ranking
# ------------------------------------------------------------
# Order by number/type of external axes for readability, then by
# internal evidence dimensions. No arithmetic priority score exists.

setorder(
  master,
  -ExternalDatabaseOverlapCount,
  -IndependentValidatedEvidenceCount,
  -CPICLevelABSupport,
  -ChEMBLSingleProteinMoASupport,
  -DGIdb_Support,
  -Tier1Metabolites,
  -UniqueMetabolites,
  -UniqueSNPs,
  Gene,
  NatMed_DrugCategory
)

# Preferred publication-facing column order.
preferred <- c(
  "Gene",
  "NatMed_DrugCategory",
  "NatMed_CategoryCode",
  "CardiovascularATCGroup",

  # Internal evidence
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
  "OverlappingGeneSNPs",
  "NearestGeneSNPs",
  "MinMRFDR",
  "MaxPPH4",
  "MinM2P",
  "SNPs",
  "Metabolites",
  "MultiMetaboliteInternalSupport",
  "MultiSNPInternalSupport",
  "FunctionalVariantInternalSupport",

  # External support overview
  "ExternalDatabaseOverlapCount",
  "IndependentValidatedEvidenceCount",
  "ExternalDatabasePattern",
  "CPICTrackedPair",
  "CPICLevelABSupport",
  "DGIdbSupport",
  "ChEMBLMoASupport",
  "ChEMBLSingleProteinMoASupport",

  # CPIC
  "CPICTrackedPair",
  "CPIC_PrimaryAB",
  "CPIC_GuidelineLinked",
  "CPIC_Drugs",
  "CPIC_Levels",

  # DGIdb
  "DGIdb_Support",
  "DGIdb_UniqueDrugs",
  "DGIdb_Drugs",
  "DGIdb_UniqueSources",
  "DGIdb_Sources",
  "DGIdb_InteractionTypes",
  "DGIdb_MaxEvidenceScore",
  "DGIdb_MaxInteractionScore",

  # ChEMBL
  "ChEMBL_Support",
  "ChEMBL_PrimarySingleProteinSupport",
  "ChEMBL_TargetCount",
  "ChEMBL_Targets",
  "ChEMBL_TargetNames",
  "ChEMBL_TargetTypes",
  "ChEMBL_ActionTypes",
  "ChEMBL_Mechanisms",
  "ChEMBL_InputDrugs"
)

keep <- intersect(preferred, names(master))
remaining <- setdiff(names(master), keep)
setcolorder(master, c(keep, remaining))

fwrite(
  master,
  master_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# ------------------------------------------------------------
# External support pattern summary
# ------------------------------------------------------------

pattern_summary <- master[
  ,
  .(
    Modules = .N,
    UniqueGenes = uniqueN(Gene),
    UniqueDrugCategories = uniqueN(NatMed_DrugCategory),
    CardiovascularModules = sum(
      CardiovascularATCGroup %in% TRUE,
      na.rm = TRUE
    ),
    ModulesWithCodingOrSpliceSupport = sum(
      FunctionalVariantInternalSupport %in% TRUE,
      na.rm = TRUE
    )
  ),
  by = .(
    ExternalDatabaseOverlapCount,
    ExternalDatabasePattern
  )
]

setorder(
  pattern_summary,
  -ExternalDatabaseOverlapCount,
  -Modules,
  ExternalDatabasePattern
)

fwrite(
  pattern_summary,
  pattern_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# ------------------------------------------------------------
# QC
# ------------------------------------------------------------

qc <- data.table(
  Metric = c(
    "Integrated gene-drug modules",
    "Unique genes",
    "Unique NatMed drug categories",

    "Modules with no external database overlap",
    "Modules with >=1 external database overlap",
    "Modules with >=2 external database overlaps",
    "Modules overlapping all 3 external databases",

    "Modules present as CPIC-tracked pairs",
    "Modules with CPIC Level A/B support",
    "Modules with CPIC guideline-linked support",

    "Modules with DGIdb support",

    "Modules with ChEMBL MoA support",
    "Modules with ChEMBL SINGLE PROTEIN MoA support",

    "Modules with CPIC + DGIdb support",
    "Modules with CPIC + ChEMBL support",
    "Modules with DGIdb + ChEMBL support",
    "Modules with CPIC + DGIdb + ChEMBL support",

    "Externally supported modules with >=2 metabolites",
    "Externally supported modules with >=2 SNPs",
    "Externally supported modules with coding/splice support"
  ),
  Value = c(
    nrow(master),
    uniqueN(master$Gene),
    uniqueN(master$NatMed_DrugCategory),

    master[ExternalDatabaseOverlapCount == 0, .N],
    master[ExternalDatabaseOverlapCount >= 1, .N],
    master[ExternalDatabaseOverlapCount >= 2, .N],
    master[ExternalDatabaseOverlapCount == 3, .N],

    master[CPICTrackedPair == TRUE, .N],
    master[CPIC_PrimaryAB == TRUE, .N],
    master[CPIC_GuidelineLinked == TRUE, .N],

    master[DGIdb_Support == TRUE, .N],

    master[ChEMBL_Support == TRUE, .N],
    master[
      ChEMBL_PrimarySingleProteinSupport == TRUE,
      .N
    ],

    master[
      CPICTrackedPair == TRUE &
        DGIdb_Support == TRUE,
      .N
    ],
    master[
      CPICTrackedPair == TRUE &
        ChEMBL_Support == TRUE,
      .N
    ],
    master[
      DGIdb_Support == TRUE &
        ChEMBL_Support == TRUE,
      .N
    ],
    master[
      CPICTrackedPair == TRUE &
        DGIdb_Support == TRUE &
        ChEMBL_Support == TRUE,
      .N
    ],

    master[
      ExternalDatabaseOverlapCount >= 1 &
        UniqueMetabolites >= 2,
      .N
    ],
    master[
      ExternalDatabaseOverlapCount >= 1 &
        UniqueSNPs >= 2,
      .N
    ],
    master[
      ExternalDatabaseOverlapCount >= 1 &
        FunctionalVariantInternalSupport == TRUE,
      .N
    ]
  )
)

fwrite(
  qc,
  qc_file,
  sep = "\t",
  quote = FALSE
)

# ------------------------------------------------------------
# Numeric formatting policy
# ------------------------------------------------------------
# Canonical master output retains MR-FDR and P-value columns as numeric
# values for downstream analysis. Scientific-notation display formatting
# will be applied only in publication-facing export tables.
#
# ------------------------------------------------------------
# Terminal summary
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("Integrated evidence summary completed successfully\n")
cat("============================================================\n\n")

print(qc)

cat("\nModules overlapping CPIC, DGIdb and ChEMBL:\n")
print(
  master[
    ExternalDatabaseOverlapCount == 3,
    .(
      Gene,
      NatMed_DrugCategory,
      UniqueSNPs,
      UniqueMetabolites,
      Tier1Metabolites,
      CPIC_Levels,
      DGIdb_UniqueSources,
      ChEMBL_TargetNames,
      ChEMBL_ActionTypes
    )
  ]
)

cat("\nModules with CPIC Level A/B or direct SINGLE PROTEIN MoA support:\n")
print(
  master[
    CPIC_PrimaryAB == TRUE |
      ChEMBL_PrimarySingleProteinSupport == TRUE,
    .(
      Gene,
      NatMed_DrugCategory,
      ExternalDatabasePattern,
      UniqueSNPs,
      UniqueMetabolites,
      Tier1Metabolites,
      CPIC_Levels,
      DGIdb_UniqueSources,
      ChEMBL_TargetNames,
      ChEMBL_ActionTypes
    )
  ]
)

cat(
  "\nOutputs:\n",
  "  ", master_file, "\n",
  "  ", pattern_file, "\n",
  "  ", qc_file, "\n\n",
  sep = ""
)

cat(
  "INTERPRETATION NOTE:\n",
  "ExternalDatabaseOverlapCount and IndependentValidatedEvidenceCount are descriptive\n",
  "counts only. They are NOT prioritisation scores and should not be used\n",
  "as evidence that one biological module is causally stronger than another.\n"
)
