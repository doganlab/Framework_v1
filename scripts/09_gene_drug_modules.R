#!/usr/bin/env Rscript

# ============================================================
# Framework_v1
# 09_gene_drug_modules.R
#
# Purpose
# -------
# Collapse the final variant-gene-metabolite-drug evidence table
# into interpretable SNP-gene-drug and gene-drug modules.
#
# IMPORTANT
# ---------
# No composite prioritisation score is calculated.
# Counts are descriptive evidence dimensions, not independent
# replications, because Nightingale metabolite traits can be
# strongly correlated.
# ============================================================

options(stringsAsFactors = FALSE, warn = 1, scipen = 999)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
source(file.path(dirname(script_path), "..", "R", "bootstrap.R"))

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("The data.table package is required.", call. = FALSE)
}
library(data.table)

input_file <- file.path(
  PROJECT_ROOT, "results", "evidence",
  "VariantGeneMetaboliteDrugEvidence.tsv"
)

output_dir <- file.path(
  PROJECT_ROOT, "results", "gene_drug_modules"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

snp_gene_drug_file <- file.path(
  output_dir, "SNP_Gene_DrugModules.tsv"
)
gene_drug_file <- file.path(
  output_dir, "Gene_DrugModules.tsv"
)
cv_gene_drug_file <- file.path(
  output_dir, "Cardiovascular_Gene_DrugModules.tsv"
)
coding_file <- file.path(
  output_dir, "CodingSplice_GeneDrugEvidence.tsv"
)
drug_summary_file <- file.path(
  output_dir, "DrugCategory_ModuleSummary.tsv"
)
gene_summary_file <- file.path(
  output_dir, "Gene_ModuleSummary.tsv"
)
qc_file <- file.path(
  output_dir, "GeneDrugModules_QC.tsv"
)

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file, call. = FALSE)
}

x <- fread(
  input_file,
  na.strings = c("", "NA"),
  showProgress = TRUE
)

required <- c(
  "SNP", "Metabolite", "DrugCategory", "CategoryCode",
  "AssignedGene", "GeneAssignmentType", "EvidenceTier",
  "MR_beta", "MR_FDR", "PP_H4", "M2_Beta", "M2_P",
  "CardiovascularATCGroup", "CodingVariant", "SpliceRelated",
  "Consequence", "IMPACT"
)

missing <- setdiff(required, names(x))
if (length(missing) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing, collapse = ", "),
    call. = FALSE
  )
}

x <- x[
  !is.na(AssignedGene) &
    AssignedGene != "" &
    EvidenceTier %in% c("Tier1_Robust", "Tier2_HighConfidence")
]

# Stable row-level flags.
x[, Tier1Flag := EvidenceTier == "Tier1_Robust"]
x[, Tier2Flag := EvidenceTier == "Tier2_HighConfidence"]
x[, OverlapFlag := GeneAssignmentType == "overlapping_gene"]
x[, NearestFlag := GeneAssignmentType == "nearest_gene"]
x[, CodingFlag := CodingVariant %in% TRUE]
x[, SpliceFlag := SpliceRelated %in% TRUE]
x[, HighImpactFlag := IMPACT == "HIGH"]
x[, ModerateImpactFlag := IMPACT == "MODERATE"]

# Helper for safe collapse of character vectors.
collapse_unique <- function(v, max_n = Inf) {
  v <- sort(unique(v[!is.na(v) & v != ""]))
  if (length(v) == 0) return(NA_character_)
  if (length(v) > max_n) {
    return(paste0(
      paste(v[seq_len(max_n)], collapse = ";"),
      ";...(+", length(v) - max_n, ")"
    ))
  }
  paste(v, collapse = ";")
}

# ------------------------------------------------------------
# 1. SNP-gene-drug modules
# ------------------------------------------------------------
# A single SNP-gene-drug row can be supported by multiple
# metabolite traits. We retain metabolite diversity explicitly.

snp_gene_drug <- x[
  ,
  .(
    UniqueMetabolites = uniqueN(Metabolite),
    Tier1Metabolites = uniqueN(Metabolite[Tier1Flag]),
    Tier2Metabolites = uniqueN(Metabolite[Tier2Flag]),
    Metabolites = collapse_unique(Metabolite, 40),
    MinMRFDR = min(MR_FDR, na.rm = TRUE),
    MaxPPH4 = max(PP_H4, na.rm = TRUE),
    MinM2P = min(M2_P, na.rm = TRUE),
    AnyCoding = any(CodingFlag, na.rm = TRUE),
    AnySpliceRelated = any(SpliceFlag, na.rm = TRUE),
    AnyHighImpact = any(HighImpactFlag, na.rm = TRUE),
    AnyModerateImpact = any(ModerateImpactFlag, na.rm = TRUE),
    AnyOverlappingGene = any(OverlapFlag, na.rm = TRUE),
    AnyNearestGene = any(NearestFlag, na.rm = TRUE),
    Consequences = collapse_unique(Consequence, 20),
    Impacts = collapse_unique(IMPACT, 10),
    CardiovascularATCGroup = any(
      CardiovascularATCGroup %in% TRUE, na.rm = TRUE
    )
  ),
  by = .(
    SNP,
    AssignedGene,
    DrugCategory,
    CategoryCode
  )
]

setorder(
  snp_gene_drug,
  -Tier1Metabolites,
  -UniqueMetabolites,
  -AnyCoding,
  -AnySpliceRelated,
  AssignedGene,
  DrugCategory,
  SNP
)

fwrite(
  snp_gene_drug,
  snp_gene_drug_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# ------------------------------------------------------------
# 2. Gene-drug modules
# ------------------------------------------------------------
# We count unique SNPs and unique metabolites separately.
# Raw chain count is deliberately NOT used as a strength score.

gene_drug <- x[
  ,
  .(
    UniqueSNPs = uniqueN(SNP),
    UniqueMetabolites = uniqueN(Metabolite),
    Tier1SNPs = uniqueN(SNP[Tier1Flag]),
    Tier1Metabolites = uniqueN(Metabolite[Tier1Flag]),
    Tier2SNPs = uniqueN(SNP[Tier2Flag]),
    Tier2Metabolites = uniqueN(Metabolite[Tier2Flag]),
    CodingSNPs = uniqueN(SNP[CodingFlag]),
    SpliceRelatedSNPs = uniqueN(SNP[SpliceFlag]),
    HighImpactSNPs = uniqueN(SNP[HighImpactFlag]),
    ModerateImpactSNPs = uniqueN(SNP[ModerateImpactFlag]),
    OverlappingGeneSNPs = uniqueN(SNP[OverlapFlag]),
    NearestGeneSNPs = uniqueN(SNP[NearestFlag]),
    MinMRFDR = min(MR_FDR, na.rm = TRUE),
    MaxPPH4 = max(PP_H4, na.rm = TRUE),
    MinM2P = min(M2_P, na.rm = TRUE),
    SNPs = collapse_unique(SNP, 50),
    Metabolites = collapse_unique(Metabolite, 50),
    Consequences = collapse_unique(Consequence, 30),
    CardiovascularATCGroup = any(
      CardiovascularATCGroup %in% TRUE, na.rm = TRUE
    )
  ),
  by = .(
    AssignedGene,
    DrugCategory,
    CategoryCode
  )
]

# Descriptive flags only; these are not scores.
gene_drug[, MultiMetaboliteModule := UniqueMetabolites >= 2]
gene_drug[, MultiSNPModule := UniqueSNPs >= 2]
gene_drug[, FunctionalVariantSupport :=
  CodingSNPs > 0 | SpliceRelatedSNPs > 0
]

setorder(
  gene_drug,
  -Tier1Metabolites,
  -UniqueMetabolites,
  -Tier1SNPs,
  -UniqueSNPs,
  -FunctionalVariantSupport,
  AssignedGene,
  DrugCategory
)

fwrite(
  gene_drug,
  gene_drug_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

fwrite(
  gene_drug[CardiovascularATCGroup %in% TRUE],
  cv_gene_drug_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# ------------------------------------------------------------
# 3. Functional variant subset
# ------------------------------------------------------------

coding_splice <- x[
  CodingFlag | SpliceFlag | HighImpactFlag | ModerateImpactFlag,
  .(
    SNP,
    AssignedGene,
    GeneAssignmentType,
    Metabolite,
    DrugCategory,
    CategoryCode,
    EvidenceTier,
    Consequence,
    IMPACT,
    CodingVariant,
    SpliceRelated,
    MR_beta,
    MR_FDR,
    PP_H4,
    M2_Beta,
    M2_P,
    CardiovascularATCGroup
  )
]

coding_splice <- unique(coding_splice)

setorder(
  coding_splice,
  -CardiovascularATCGroup,
  EvidenceTier,
  AssignedGene,
  DrugCategory,
  Metabolite,
  SNP
)

fwrite(
  coding_splice,
  coding_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# ------------------------------------------------------------
# 4. Drug- and gene-level descriptive summaries
# ------------------------------------------------------------

drug_summary <- gene_drug[
  ,
  .(
    GeneDrugModules = .N,
    UniqueGenes = uniqueN(AssignedGene),
    GenesWithMultiMetaboliteSupport = uniqueN(
      AssignedGene[MultiMetaboliteModule == TRUE]
    ),
    GenesWithMultiSNPSupport = uniqueN(
      AssignedGene[MultiSNPModule == TRUE]
    ),
    GenesWithFunctionalVariantSupport = uniqueN(
      AssignedGene[FunctionalVariantSupport == TRUE]
    ),
    TotalUniqueGeneSNPpairs = sum(UniqueSNPs),
    MaxMetabolitesPerGene = max(UniqueMetabolites),
    MaxSNPsPerGene = max(UniqueSNPs),
    CardiovascularATCGroup = any(
      CardiovascularATCGroup %in% TRUE, na.rm = TRUE
    )
  ),
  by = .(DrugCategory, CategoryCode)
]

setorder(
  drug_summary,
  -CardiovascularATCGroup,
  -GenesWithMultiMetaboliteSupport,
  -UniqueGenes,
  DrugCategory
)

fwrite(
  drug_summary,
  drug_summary_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

gene_summary <- gene_drug[
  ,
  .(
    DrugCategories = uniqueN(DrugCategory),
    CardiovascularDrugCategories = uniqueN(
      DrugCategory[CardiovascularATCGroup %in% TRUE]
    ),
    TotalUniqueSNPDrugPairs = sum(UniqueSNPs),
    MaxMetabolitesInOneDrugModule = max(UniqueMetabolites),
    ModulesWithMultiMetaboliteSupport = sum(MultiMetaboliteModule),
    ModulesWithFunctionalVariantSupport = sum(FunctionalVariantSupport)
  ),
  by = AssignedGene
]

setorder(
  gene_summary,
  -CardiovascularDrugCategories,
  -ModulesWithMultiMetaboliteSupport,
  -DrugCategories,
  AssignedGene
)

fwrite(
  gene_summary,
  gene_summary_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# ------------------------------------------------------------
# 5. QC
# ------------------------------------------------------------

qc <- data.table(
  Metric = c(
    "Input evidence-chain rows",
    "Unique SNPs",
    "Unique assigned genes",
    "Unique metabolites",
    "Unique drug categories",
    "SNP-gene-drug modules",
    "Gene-drug modules",
    "Cardiovascular gene-drug modules",
    "Gene-drug modules with >=2 metabolites",
    "Gene-drug modules with >=2 SNPs",
    "Gene-drug modules with coding/splice support",
    "Gene-drug modules with HIGH-impact SNP support",
    "Unique coding/splice evidence rows"
  ),
  Value = c(
    nrow(x),
    uniqueN(x$SNP),
    uniqueN(x$AssignedGene),
    uniqueN(x$Metabolite),
    uniqueN(x$DrugCategory),
    nrow(snp_gene_drug),
    nrow(gene_drug),
    gene_drug[CardiovascularATCGroup %in% TRUE, .N],
    gene_drug[MultiMetaboliteModule == TRUE, .N],
    gene_drug[MultiSNPModule == TRUE, .N],
    gene_drug[FunctionalVariantSupport == TRUE, .N],
    gene_drug[HighImpactSNPs > 0, .N],
    nrow(coding_splice)
  )
)

fwrite(
  qc,
  qc_file,
  sep = "\t",
  quote = FALSE
)

cat("\n============================================================\n")
cat("Gene-drug module analysis completed successfully\n")
cat("============================================================\n\n")
print(qc)

cat(
  "\nMain outputs:\n",
  "  ", gene_drug_file, "\n",
  "  ", snp_gene_drug_file, "\n",
  "  ", cv_gene_drug_file, "\n",
  "  ", coding_file, "\n",
  "\nSummaries:\n",
  "  ", drug_summary_file, "\n",
  "  ", gene_summary_file, "\n",
  "  ", qc_file, "\n",
  sep = ""
)

cat(
  "\nNOTE: module counts are descriptive evidence dimensions. ",
  "They must not be interpreted as independent replications or ",
  "as a composite prioritisation score.\n"
)
