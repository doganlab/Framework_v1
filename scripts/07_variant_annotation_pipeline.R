#!/usr/bin/env Rscript

options(
  stringsAsFactors = FALSE,
  scipen = 999
)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
source(file.path(dirname(script_path), "..", "R", "bootstrap.R"))
setwd(PROJECT_ROOT)

# ============================================================================
# Framework_v1
# 07_variant_annotation_pipeline.R
#
# Purpose:
#   Validate parsed VEP annotations, generate QC summaries and prepare the
#   final variant-gene evidence table for downstream drug and PGx integration.
#
# Run from Framework_v1 root:
#   Rscript scripts/07_variant_annotation_pipeline.R
# ============================================================================


# ----------------------------------------------------------------------------
# 1. Root-directory check
# ----------------------------------------------------------------------------

required_root_dirs <- c(
  "scripts",
  "data",
  "results"
)

missing_root_dirs <- required_root_dirs[
  !dir.exists(required_root_dirs)
]

if (length(missing_root_dirs) > 0) {
  stop(
    paste0(
      "This script must be run from the Framework_v1 root directory.\n",
      "Missing directories: ",
      paste(missing_root_dirs, collapse = ", "),
      "\nCurrent working directory: ",
      getwd()
    )
  )
}

cat("\nFramework root:\n  ", getwd(), "\n\n", sep = "")


# ----------------------------------------------------------------------------
# 2. Input and output paths
# ----------------------------------------------------------------------------

annotation_dir <- "results/variant_annotation"

association_file <- file.path(
  annotation_dir,
  "VariantGeneAssociations.tsv"
)

summary_file <- file.path(
  annotation_dir,
  "VariantGeneSummary.tsv"
)

merged_file <- file.path(
  annotation_dir,
  "InstrumentVariantGeneEvidence.tsv"
)

final_file <- file.path(
  annotation_dir,
  "VariantGeneEvidence.tsv"
)

qc_dir <- file.path(
  annotation_dir,
  "qc"
)

dir.create(
  qc_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

qc_summary_file <- file.path(
  qc_dir,
  "VariantAnnotation_QC_Summary.tsv"
)

impact_summary_file <- file.path(
  qc_dir,
  "ImpactSummary.tsv"
)

consequence_summary_file <- file.path(
  qc_dir,
  "ConsequenceSummary.tsv"
)

gene_assignment_summary_file <- file.path(
  qc_dir,
  "GeneAssignmentSummary.tsv"
)

metabolite_summary_file <- file.path(
  qc_dir,
  "MetaboliteVariantGeneSummary.tsv"
)

unannotated_file <- file.path(
  qc_dir,
  "UnannotatedVariants.tsv"
)

multi_gene_file <- file.path(
  qc_dir,
  "MultiGeneVariants.tsv"
)

duplicate_file <- file.path(
  qc_dir,
  "DuplicateInstrumentRows.tsv"
)

pipeline_report_file <- file.path(
  qc_dir,
  "VariantAnnotation_PipelineReport.txt"
)


# ----------------------------------------------------------------------------
# 3. Helper functions
# ----------------------------------------------------------------------------

check_file <- function(path) {
  
  if (!file.exists(path)) {
    stop("Required input file not found: ", path)
  }
  
  file_size <- file.info(path)$size
  
  if (is.na(file_size) || file_size == 0) {
    stop("Input file is empty: ", path)
  }
}


read_tsv <- function(path) {
  
  read.delim(
    path,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    check.names = FALSE,
    na.strings = c("", "NA"),
    stringsAsFactors = FALSE
  )
}


write_tsv <- function(data, path) {
  
  write.table(
    data,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = ""
  )
}


safe_unique_count <- function(x) {
  
  length(
    unique(
      x[
        !is.na(x) &
          nzchar(trimws(as.character(x)))
      ]
    )
  )
}


normalise_character <- function(x) {
  
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  x
}


make_frequency_table <- function(x, column_name) {
  
  x <- normalise_character(x)
  x[is.na(x)] <- "Unassigned"
  
  result <- as.data.frame(
    table(x),
    stringsAsFactors = FALSE
  )
  
  colnames(result) <- c(
    column_name,
    "Count"
  )
  
  result$Percent <- round(
    100 * result$Count / sum(result$Count),
    2
  )
  
  result <- result[
    order(
      -result$Count,
      result[[column_name]]
    ),
    ,
    drop = FALSE
  ]
  
  rownames(result) <- NULL
  
  result
}


# ----------------------------------------------------------------------------
# 4. Validate input files
# ----------------------------------------------------------------------------

input_files <- c(
  association_file,
  summary_file,
  merged_file
)

invisible(
  lapply(
    input_files,
    check_file
  )
)

cat("Input files found:\n")

for (path in input_files) {
  cat("  ", path, "\n", sep = "")
}

cat("\n")


# ----------------------------------------------------------------------------
# 5. Read data
# ----------------------------------------------------------------------------

associations <- read_tsv(
  association_file
)

variant_summary <- read_tsv(
  summary_file
)

evidence <- read_tsv(
  merged_file
)

cat("Rows loaded:\n")
cat("  VariantGeneAssociations: ", nrow(associations), "\n", sep = "")
cat("  VariantGeneSummary: ", nrow(variant_summary), "\n", sep = "")
cat("  InstrumentVariantGeneEvidence: ", nrow(evidence), "\n\n", sep = "")


# ----------------------------------------------------------------------------
# 6. Required-column checks
# ----------------------------------------------------------------------------

required_association_columns <- c(
  "SNP",
  "Consequence",
  "IMPACT"
)

required_summary_columns <- c(
  "SNP",
  "AssignedGene",
  "GeneAssignmentType",
  "Consequence",
  "IMPACT"
)

required_evidence_columns <- c(
  "Metabolite",
  "StudyAccession",
  "SNP",
  "CHR",
  "BP",
  "EffectAllele",
  "OtherAllele",
  "BetaExposure",
  "SEExposure",
  "PExposure",
  "BetaOutcome",
  "SEOutcome",
  "POutcome",
  "FStatistic",
  "AssignedGene",
  "GeneAssignmentType",
  "Consequence",
  "IMPACT"
)

# R2 is intentionally NOT required in the current NatMed pipeline.
# The present InstrumentEvidence table carries per-SNP F statistics,
# while SNP-level exposure R2 was not estimated because the required
# trait-specific sample-size information is not propagated uniformly
# across all mGWAS traits. If an R2 column is present in a future run,
# it will still be retained automatically in the final output.

check_columns <- function(data, required, label) {
  
  missing <- setdiff(
    required,
    colnames(data)
  )
  
  if (length(missing) > 0) {
    stop(
      label,
      " is missing required columns: ",
      paste(missing, collapse = ", ")
    )
  }
}

check_columns(
  associations,
  required_association_columns,
  "VariantGeneAssociations.tsv"
)

check_columns(
  variant_summary,
  required_summary_columns,
  "VariantGeneSummary.tsv"
)

check_columns(
  evidence,
  required_evidence_columns,
  "InstrumentVariantGeneEvidence.tsv"
)


# ----------------------------------------------------------------------------
# 7. Normalise important fields
# ----------------------------------------------------------------------------

character_columns <- intersect(
  c(
    "Metabolite",
    "StudyAccession",
    "SNP",
    "EffectAllele",
    "OtherAllele",
    "AssignedGene",
    "GeneAssignmentType",
    "Gene",
    "Feature",
    "Feature_type",
    "Consequence",
    "IMPACT",
    "BIOTYPE",
    "CANONICAL",
    "VARIANT_CLASS",
    "AllAssignedGenes",
    "NEAREST"
  ),
  colnames(evidence)
)

for (column in character_columns) {
  evidence[[column]] <- normalise_character(
    evidence[[column]]
  )
}

evidence$IMPACT <- toupper(
  evidence$IMPACT
)

evidence$EffectAllele <- toupper(
  evidence$EffectAllele
)

evidence$OtherAllele <- toupper(
  evidence$OtherAllele
)


# ----------------------------------------------------------------------------
# 8. Annotation status
# ----------------------------------------------------------------------------

evidence$VEPAnnotated <- (
  !is.na(evidence$Consequence) |
    !is.na(evidence$AssignedGene)
)

evidence$GeneAssigned <- (
  !is.na(evidence$AssignedGene) &
    nzchar(evidence$AssignedGene)
)

evidence$CanonicalTranscript <- (
  !is.na(evidence$CANONICAL) &
    evidence$CANONICAL == "YES"
)

evidence$ProteinCoding <- (
  !is.na(evidence$BIOTYPE) &
    evidence$BIOTYPE == "protein_coding"
)

evidence$CodingVariant <- grepl(
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
  ifelse(
    is.na(evidence$Consequence),
    "",
    evidence$Consequence
  )
)

evidence$SpliceRelated <- grepl(
  "splice",
  ifelse(
    is.na(evidence$Consequence),
    "",
    evidence$Consequence
  ),
  ignore.case = TRUE
)

evidence$RegulatoryProximity <- grepl(
  paste(
    c(
      "upstream_gene_variant",
      "downstream_gene_variant",
      "5_prime_UTR_variant",
      "3_prime_UTR_variant"
    ),
    collapse = "|"
  ),
  ifelse(
    is.na(evidence$Consequence),
    "",
    evidence$Consequence
  )
)


# ----------------------------------------------------------------------------
# 9. Functional annotation class
# ----------------------------------------------------------------------------

evidence$FunctionalClass <- "Other"

evidence$FunctionalClass[
  grepl(
    "intergenic_variant",
    ifelse(
      is.na(evidence$Consequence),
      "",
      evidence$Consequence
    )
  )
] <- "Intergenic"

evidence$FunctionalClass[
  grepl(
    "upstream_gene_variant|downstream_gene_variant",
    ifelse(
      is.na(evidence$Consequence),
      "",
      evidence$Consequence
    )
  )
] <- "Gene proximal"

evidence$FunctionalClass[
  grepl(
    "intron_variant",
    ifelse(
      is.na(evidence$Consequence),
      "",
      evidence$Consequence
    )
  )
] <- "Intronic"

evidence$FunctionalClass[
  grepl(
    "UTR_variant",
    ifelse(
      is.na(evidence$Consequence),
      "",
      evidence$Consequence
    )
  )
] <- "UTR"

evidence$FunctionalClass[
  grepl(
    "non_coding_transcript",
    ifelse(
      is.na(evidence$Consequence),
      "",
      evidence$Consequence
    )
  )
] <- "Non-coding transcript"

evidence$FunctionalClass[
  evidence$SpliceRelated
] <- "Splice related"

evidence$FunctionalClass[
  evidence$CodingVariant
] <- "Coding"


# ----------------------------------------------------------------------------
# 10. Annotation priority score
#
# This is not a biological causal score. It is only an annotation-priority
# index for sorting variants within the Framework.
# ----------------------------------------------------------------------------

impact_score <- c(
  "HIGH" = 4,
  "MODERATE" = 3,
  "LOW" = 2,
  "MODIFIER" = 1
)

evidence$ImpactScore <- unname(
  impact_score[evidence$IMPACT]
)

evidence$ImpactScore[
  is.na(evidence$ImpactScore)
] <- 0

evidence$AnnotationPriorityScore <- (
  evidence$ImpactScore +
    ifelse(evidence$GeneAssigned, 1, 0) +
    ifelse(evidence$CanonicalTranscript, 1, 0) +
    ifelse(evidence$ProteinCoding, 1, 0) +
    ifelse(evidence$CodingVariant, 2, 0) +
    ifelse(evidence$SpliceRelated, 2, 0)
)


# ----------------------------------------------------------------------------
# 11. Variant-level QC
# ----------------------------------------------------------------------------

unique_instrument_snps <- unique(
  evidence$SNP
)

unique_summary_snps <- unique(
  variant_summary$SNP
)

unannotated_snps <- setdiff(
  unique_instrument_snps,
  unique_summary_snps
)

unannotated <- evidence[
  evidence$SNP %in% unannotated_snps |
    !evidence$VEPAnnotated,
  ,
  drop = FALSE
]

write_tsv(
  unannotated,
  unannotated_file
)


# ----------------------------------------------------------------------------
# 12. Multi-gene variants
# ----------------------------------------------------------------------------

if ("AssignedGeneCount" %in% colnames(variant_summary)) {
  
  variant_summary$AssignedGeneCount <- suppressWarnings(
    as.integer(variant_summary$AssignedGeneCount)
  )
  
  multi_gene <- variant_summary[
    !is.na(variant_summary$AssignedGeneCount) &
      variant_summary$AssignedGeneCount > 1,
    ,
    drop = FALSE
  ]
  
} else if ("AllAssignedGenes" %in% colnames(variant_summary)) {
  
  gene_counts <- vapply(
    strsplit(
      ifelse(
        is.na(variant_summary$AllAssignedGenes),
        "",
        variant_summary$AllAssignedGenes
      ),
      ";",
      fixed = TRUE
    ),
    function(x) {
      x <- x[nzchar(x)]
      length(unique(x))
    },
    integer(1)
  )
  
  variant_summary$AssignedGeneCount <- gene_counts
  
  multi_gene <- variant_summary[
    gene_counts > 1,
    ,
    drop = FALSE
  ]
  
} else {
  
  multi_gene <- variant_summary[0, , drop = FALSE]
}

write_tsv(
  multi_gene,
  multi_gene_file
)


# ----------------------------------------------------------------------------
# 13. Duplicate instrument-row check
# ----------------------------------------------------------------------------

duplicate_key <- paste(
  evidence$Metabolite,
  evidence$StudyAccession,
  evidence$SNP,
  sep = "|"
)

duplicate_rows <- evidence[
  duplicated(duplicate_key) |
    duplicated(duplicate_key, fromLast = TRUE),
  ,
  drop = FALSE
]

write_tsv(
  duplicate_rows,
  duplicate_file
)


# ----------------------------------------------------------------------------
# 14. Frequency summaries
# ----------------------------------------------------------------------------

impact_summary <- make_frequency_table(
  evidence$IMPACT,
  "IMPACT"
)

write_tsv(
  impact_summary,
  impact_summary_file
)


consequence_summary <- make_frequency_table(
  evidence$Consequence,
  "Consequence"
)

write_tsv(
  consequence_summary,
  consequence_summary_file
)


gene_assignment_summary <- make_frequency_table(
  evidence$GeneAssignmentType,
  "GeneAssignmentType"
)

write_tsv(
  gene_assignment_summary,
  gene_assignment_summary_file
)


# ----------------------------------------------------------------------------
# 15. Metabolite-level summary
# ----------------------------------------------------------------------------

metabolites <- sort(
  unique(
    evidence$Metabolite[
      !is.na(evidence$Metabolite)
    ]
  )
)

metabolite_summary <- do.call(
  rbind,
  lapply(
    metabolites,
    function(metabolite_name) {
      
      x <- evidence[
        evidence$Metabolite == metabolite_name,
        ,
        drop = FALSE
      ]
      
      data.frame(
        Metabolite = metabolite_name,
        StudyAccession = paste(
          sort(
            unique(
              x$StudyAccession[
                !is.na(x$StudyAccession)
              ]
            )
          ),
          collapse = ";"
        ),
        InstrumentRows = nrow(x),
        UniqueSNPs = safe_unique_count(x$SNP),
        AnnotatedSNPs = safe_unique_count(
          x$SNP[x$VEPAnnotated]
        ),
        GeneAssignedSNPs = safe_unique_count(
          x$SNP[x$GeneAssigned]
        ),
        UniqueAssignedGenes = safe_unique_count(
          x$AssignedGene
        ),
        CodingSNPs = safe_unique_count(
          x$SNP[x$CodingVariant]
        ),
        SpliceRelatedSNPs = safe_unique_count(
          x$SNP[x$SpliceRelated]
        ),
        CanonicalTranscriptSNPs = safe_unique_count(
          x$SNP[x$CanonicalTranscript]
        ),
        HighImpactSNPs = safe_unique_count(
          x$SNP[x$IMPACT == "HIGH"]
        ),
        ModerateImpactSNPs = safe_unique_count(
          x$SNP[x$IMPACT == "MODERATE"]
        ),
        LowImpactSNPs = safe_unique_count(
          x$SNP[x$IMPACT == "LOW"]
        ),
        ModifierSNPs = safe_unique_count(
          x$SNP[x$IMPACT == "MODIFIER"]
        ),
        stringsAsFactors = FALSE
      )
    }
  )
)

rownames(metabolite_summary) <- NULL

write_tsv(
  metabolite_summary,
  metabolite_summary_file
)


# ----------------------------------------------------------------------------
# 16. Select and order final evidence columns
# ----------------------------------------------------------------------------

preferred_columns <- c(
  "Metabolite",
  "StudyAccession",
  "SNP",
  "CHR",
  "BP",
  "EffectAllele",
  "OtherAllele",
  "EffectAlleleFrequency",
  "BetaExposure",
  "SEExposure",
  "PExposure",
  "BetaOutcome",
  "SEOutcome",
  "POutcome",
  "FStatistic",
  "R2",
  "InstrumentStrengthClass",
  "WeakInstrument",
  "HarmonizationStatus",
  "Palindromic",
  "Ambiguous",
  "MHC",
  "SteigerDirection",
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
  "TranscriptCount",
  "VEPAnnotated",
  "GeneAssigned",
  "CanonicalTranscript",
  "ProteinCoding",
  "CodingVariant",
  "SpliceRelated",
  "RegulatoryProximity",
  "FunctionalClass",
  "ImpactScore",
  "AnnotationPriorityScore",
  "SourceFile"
)

final_columns <- c(
  intersect(
    preferred_columns,
    colnames(evidence)
  ),
  setdiff(
    colnames(evidence),
    preferred_columns
  )
)

final_evidence <- evidence[
  ,
  final_columns,
  drop = FALSE
]

final_evidence <- final_evidence[
  order(
    final_evidence$Metabolite,
    -final_evidence$AnnotationPriorityScore,
    final_evidence$PExposure,
    final_evidence$SNP,
    na.last = TRUE
  ),
  ,
  drop = FALSE
]

rownames(final_evidence) <- NULL

write_tsv(
  final_evidence,
  final_file
)


# ----------------------------------------------------------------------------
# 17. QC summary
# ----------------------------------------------------------------------------

qc_summary <- data.frame(
  Metric = c(
    "Instrument evidence rows",
    "Unique instrument SNPs",
    "VEP summary SNPs",
    "VEP association rows",
    "Annotated instrument SNPs",
    "Unannotated instrument SNPs",
    "Instrument SNPs with assigned gene",
    "Unique assigned genes",
    "Multi-gene SNPs",
    "Canonical transcript SNPs",
    "Protein-coding SNPs",
    "Coding SNPs",
    "Splice-related SNPs",
    "HIGH impact SNPs",
    "MODERATE impact SNPs",
    "LOW impact SNPs",
    "MODIFIER impact SNPs",
    "Duplicate metabolite-study-SNP rows",
    "Metabolites"
  ),
  Value = c(
    nrow(evidence),
    safe_unique_count(evidence$SNP),
    safe_unique_count(variant_summary$SNP),
    nrow(associations),
    safe_unique_count(
      evidence$SNP[evidence$VEPAnnotated]
    ),
    safe_unique_count(
      evidence$SNP[!evidence$VEPAnnotated]
    ),
    safe_unique_count(
      evidence$SNP[evidence$GeneAssigned]
    ),
    safe_unique_count(evidence$AssignedGene),
    safe_unique_count(multi_gene$SNP),
    safe_unique_count(
      evidence$SNP[evidence$CanonicalTranscript]
    ),
    safe_unique_count(
      evidence$SNP[evidence$ProteinCoding]
    ),
    safe_unique_count(
      evidence$SNP[evidence$CodingVariant]
    ),
    safe_unique_count(
      evidence$SNP[evidence$SpliceRelated]
    ),
    safe_unique_count(
      evidence$SNP[evidence$IMPACT == "HIGH"]
    ),
    safe_unique_count(
      evidence$SNP[evidence$IMPACT == "MODERATE"]
    ),
    safe_unique_count(
      evidence$SNP[evidence$IMPACT == "LOW"]
    ),
    safe_unique_count(
      evidence$SNP[evidence$IMPACT == "MODIFIER"]
    ),
    nrow(duplicate_rows),
    safe_unique_count(evidence$Metabolite)
  ),
  stringsAsFactors = FALSE
)

write_tsv(
  qc_summary,
  qc_summary_file
)


# ----------------------------------------------------------------------------
# 18. Text report
# ----------------------------------------------------------------------------

report_lines <- c(
  "Framework_v1 Variant Annotation Pipeline",
  "========================================",
  "",
  paste("Run date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste("Working directory:", getwd()),
  "",
  "Input files:",
  paste0("  ", association_file),
  paste0("  ", summary_file),
  paste0("  ", merged_file),
  "",
  "Main results:",
  paste0(
    "  Instrument rows: ",
    nrow(evidence)
  ),
  paste0(
    "  Unique instrument SNPs: ",
    safe_unique_count(evidence$SNP)
  ),
  paste0(
    "  Annotated SNPs: ",
    safe_unique_count(
      evidence$SNP[evidence$VEPAnnotated]
    )
  ),
  paste0(
    "  Gene-assigned SNPs: ",
    safe_unique_count(
      evidence$SNP[evidence$GeneAssigned]
    )
  ),
  paste0(
    "  Unique assigned genes: ",
    safe_unique_count(evidence$AssignedGene)
  ),
  paste0(
    "  Coding SNPs: ",
    safe_unique_count(
      evidence$SNP[evidence$CodingVariant]
    )
  ),
  paste0(
    "  Splice-related SNPs: ",
    safe_unique_count(
      evidence$SNP[evidence$SpliceRelated]
    )
  ),
  paste0(
    "  Multi-gene SNPs: ",
    safe_unique_count(multi_gene$SNP)
  ),
  paste0(
    "  Unannotated SNPs: ",
    safe_unique_count(unannotated$SNP)
  ),
  "",
  "Important note:",
  paste0(
    "  AnnotationPriorityScore is an internal ranking index and must not ",
    "be interpreted as evidence of causality or clinical actionability."
  ),
  "",
  "Outputs:",
  paste0("  ", final_file),
  paste0("  ", qc_summary_file),
  paste0("  ", impact_summary_file),
  paste0("  ", consequence_summary_file),
  paste0("  ", gene_assignment_summary_file),
  paste0("  ", metabolite_summary_file),
  paste0("  ", unannotated_file),
  paste0("  ", multi_gene_file),
  paste0("  ", duplicate_file)
)

writeLines(
  report_lines,
  con = pipeline_report_file
)


# ----------------------------------------------------------------------------
# 19. Console report
# ----------------------------------------------------------------------------

cat("\n============================================================\n")
cat("Variant annotation pipeline completed successfully\n")
cat("============================================================\n\n")

print(
  qc_summary,
  row.names = FALSE
)

cat("\nMain output:\n")
cat("  ", final_file, "\n", sep = "")

cat("\nQC outputs:\n")
cat("  ", qc_summary_file, "\n", sep = "")
cat("  ", metabolite_summary_file, "\n", sep = "")
cat("  ", pipeline_report_file, "\n", sep = "")

cat("\n")
