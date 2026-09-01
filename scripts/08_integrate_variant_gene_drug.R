#!/usr/bin/env Rscript

# ============================================================
# Framework_v1
# 08_integrate_variant_gene_drug.R
#
# Nature Medicine version
#
# Purpose
# -------
# Add VEP-derived variant-gene annotations to the existing
# SNP -> metabolite -> drug-category evidence chains WITHOUT
# changing their row structure.
#
# Input:
#   data/processed/VariantDrugClassEvidence.tsv
#   results/variant_annotation/VariantGeneSummary.tsv
#
# Output:
#   results/evidence/VariantGeneMetaboliteDrugEvidence.tsv
# ============================================================

options(stringsAsFactors = FALSE, warn = 1, scipen = 999)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
source(file.path(dirname(script_path), "..", "R", "bootstrap.R"))

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("The data.table package is required.", call. = FALSE)
}
library(data.table)

variant_drug_file <- VARIANT_DRUG_EVIDENCE_FILE

variant_gene_file <- file.path(
  PROJECT_ROOT,
  "results",
  "variant_annotation",
  "VariantGeneSummary.tsv"
)

output_dir <- file.path(
  PROJECT_ROOT,
  "results",
  "evidence"
)

output_file <- file.path(
  output_dir,
  "VariantGeneMetaboliteDrugEvidence.tsv"
)

qc_file <- file.path(
  output_dir,
  "VariantGeneMetaboliteDrugIntegration_QC.tsv"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(
  variant_drug_file,
  variant_gene_file
)

missing <- required_files[!file.exists(required_files)]

if (length(missing) > 0) {
  stop(
    "Missing required file(s):\n",
    paste(missing, collapse = "\n"),
    call. = FALSE
  )
}

cat("\n=== Variant-gene-metabolite-drug integration ===\n")

drug <- fread(
  variant_drug_file,
  na.strings = c("", "NA"),
  showProgress = TRUE
)

gene <- fread(
  variant_gene_file,
  na.strings = c("", "NA"),
  showProgress = FALSE
)

if (!"SNP" %in% names(drug)) {
  stop("VariantDrugClassEvidence has no SNP column.", call. = FALSE)
}

if (!"SNP" %in% names(gene)) {
  stop("VariantGeneSummary has no SNP column.", call. = FALSE)
}

if (anyDuplicated(gene$SNP) > 0) {
  stop(
    "VariantGeneSummary contains duplicated SNP rows; expected one representative row per SNP.",
    call. = FALSE
  )
}

preferred_gene_columns <- c(
  "SNP",
  "Location",
  "Allele",
  "AssignedGene",
  "GeneAssignmentType",
  "AllAssignedGenes",
  "AssignedGeneCount",
  "Gene",
  "Feature",
  "Feature_type",
  "BIOTYPE",
  "CANONICAL",
  "Consequence",
  "IMPACT",
  "VARIANT_CLASS",
  "DISTANCE",
  "Amino_acids",
  "Codons",
  "Protein_position",
  "HGNC_ID",
  "NEAREST",
  "TranscriptCount"
)

gene_columns <- intersect(
  preferred_gene_columns,
  names(gene)
)

gene <- gene[, ..gene_columns]

# Derived annotation flags. These are descriptive/QC features,
# not evidence of drug efficacy or causal pharmacogenomics.
gene[, CanonicalTranscript :=
  !is.na(CANONICAL) & CANONICAL == "YES"
]

gene[, ProteinCoding :=
  !is.na(BIOTYPE) & BIOTYPE == "protein_coding"
]

gene[, CodingVariant := grepl(
  paste(
    c(
      "missense_variant",
      "synonymous_variant",
      "stop_gained",
      "stop_lost",
      "start_lost",
      "frameshift_variant",
      "inframe_insertion",
      "inframe_deletion",
      "protein_altering_variant",
      "coding_sequence_variant"
    ),
    collapse = "|"
  ),
  fifelse(is.na(Consequence), "", Consequence)
)]

gene[, SpliceRelated := grepl(
  "splice",
  fifelse(is.na(Consequence), "", Consequence),
  ignore.case = TRUE
)]

gene[, FunctionalClass := fifelse(
  CodingVariant,
  "Coding",
  fifelse(
    SpliceRelated,
    "Splice related",
    fifelse(
      grepl("UTR_variant", fifelse(is.na(Consequence), "", Consequence)),
      "UTR",
      fifelse(
        grepl("intron_variant", fifelse(is.na(Consequence), "", Consequence)),
        "Intronic",
        fifelse(
          grepl(
            "upstream_gene_variant|downstream_gene_variant",
            fifelse(is.na(Consequence), "", Consequence)
          ),
          "Gene proximal",
          fifelse(
            grepl(
              "intergenic_variant",
              fifelse(is.na(Consequence), "", Consequence)
            ),
            "Intergenic",
            "Other"
          )
        )
      )
    )
  )
)]

impact_score <- c(
  HIGH = 4,
  MODERATE = 3,
  LOW = 2,
  MODIFIER = 1
)

gene[, ImpactScore := unname(impact_score[IMPACT])]
gene[is.na(ImpactScore), ImpactScore := 0]

gene[, AnnotationPriorityScore := (
  ImpactScore +
    fifelse(!is.na(AssignedGene) & AssignedGene != "", 1, 0) +
    fifelse(CanonicalTranscript, 1, 0) +
    fifelse(ProteinCoding, 1, 0) +
    fifelse(CodingVariant, 2, 0) +
    fifelse(SpliceRelated, 2, 0)
)]

drug[, OriginalRowID := .I]
original_rows <- nrow(drug)

merged <- merge(
  drug,
  gene,
  by = "SNP",
  all.x = TRUE,
  sort = FALSE
)

setorder(merged, OriginalRowID)
merged[, OriginalRowID := NULL]

if (nrow(merged) != original_rows) {
  stop(
    "Row count changed after VEP merge. Before: ",
    original_rows,
    "; after: ",
    nrow(merged),
    call. = FALSE
  )
}

if (!identical(
  as.character(drug$SNP),
  as.character(merged$SNP)
)) {
  # merge() does not necessarily preserve SNP vector order even after OriginalRowID
  # because the SNP column is moved; validate using rowwise core keys instead.
  core_cols <- intersect(
    c("SNP", "Metabolite", "DrugCategory"),
    names(drug)
  )
  before_keys <- do.call(
    paste,
    c(drug[, ..core_cols], sep = "||")
  )
  after_keys <- do.call(
    paste,
    c(merged[, ..core_cols], sep = "||")
  )

  if (!identical(before_keys, after_keys)) {
    stop(
      "Original variant-metabolite-drug row order changed during integration.",
      call. = FALSE
    )
  }
}

fwrite(
  merged,
  output_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

qc <- data.table(
  Metric = c(
    "Input variant-metabolite-drug rows",
    "Output variant-gene-metabolite-drug rows",
    "Unique SNPs in evidence chains",
    "Unique SNPs in VariantGeneSummary",
    "Rows with any VEP representative annotation",
    "Rows with assigned gene",
    "Unique assigned genes",
    "Rows with overlapping-gene assignment",
    "Rows with nearest-gene assignment",
    "Rows with coding variant annotation",
    "Rows with splice-related annotation",
    "Tier1 evidence rows",
    "Tier2 evidence rows",
    "Cardiovascular-ATC evidence rows"
  ),
  Value = c(
    nrow(drug),
    nrow(merged),
    uniqueN(drug$SNP),
    uniqueN(gene$SNP),
    sum(!is.na(merged$Consequence), na.rm = TRUE),
    sum(
      !is.na(merged$AssignedGene) &
        merged$AssignedGene != "",
      na.rm = TRUE
    ),
    uniqueN(
      merged[
        !is.na(AssignedGene) & AssignedGene != "",
        AssignedGene
      ]
    ),
    sum(merged$GeneAssignmentType == "overlapping_gene", na.rm = TRUE),
    sum(merged$GeneAssignmentType == "nearest_gene", na.rm = TRUE),
    sum(merged$CodingVariant %in% TRUE, na.rm = TRUE),
    sum(merged$SpliceRelated %in% TRUE, na.rm = TRUE),
    sum(merged$EvidenceTier == "Tier1_Robust", na.rm = TRUE),
    sum(merged$EvidenceTier == "Tier2_HighConfidence", na.rm = TRUE),
    sum(merged$CardiovascularATCGroup %in% TRUE, na.rm = TRUE)
  )
)

fwrite(
  qc,
  qc_file,
  sep = "\t",
  quote = FALSE
)

cat("\nIntegration completed.\n")
print(qc)

cat(
  "\nFinal evidence: ",
  output_file,
  "\nQC: ",
  qc_file,
  "\n\n",
  sep = ""
)
