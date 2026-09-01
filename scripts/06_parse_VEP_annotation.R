#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, scipen = 999)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
source(file.path(dirname(script_path), "..", "R", "bootstrap.R"))

# -------------------------------------------------------------------------
# Paths -- current mGWAS coordinates are GRCh38
# -------------------------------------------------------------------------

vep_file <- file.path(
  DIR_INTERMEDIATE,
  "variant_annotation",
  "VEP_annotation_GRCh38.tsv"
)

output_dir <- file.path(
  PROJECT_ROOT,
  "results",
  "variant_annotation"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

long_output <- file.path(
  output_dir,
  "VariantGeneAssociations.tsv"
)

summary_output <- file.path(
  output_dir,
  "VariantGeneSummary.tsv"
)

merged_output <- file.path(
  output_dir,
  "InstrumentVariantGeneEvidence.tsv"
)

instrument_file <- INSTRUMENT_EVIDENCE_FILE

if (!file.exists(vep_file)) {
  stop("VEP output not found: ", vep_file, call. = FALSE)
}

if (!file.exists(instrument_file)) {
  stop("InstrumentEvidence.tsv not found: ", instrument_file, call. = FALSE)
}

# -------------------------------------------------------------------------
# Read VEP tabular output
# -------------------------------------------------------------------------

read_vep_table <- function(path) {
  
  lines <- readLines(
    path,
    warn = FALSE
  )
  
  header_index <- grep(
    "^#Uploaded_variation",
    lines
  )
  
  if (length(header_index) == 0) {
    stop("VEP column header was not found.")
  }
  
  header_index <- header_index[1]
  
  header <- sub(
    "^#",
    "",
    lines[header_index]
  )
  
  data_lines <- lines[
    seq.int(
      from = header_index + 1,
      to = length(lines)
    )
  ]
  
  data_lines <- data_lines[
    nzchar(data_lines) &
      !grepl("^#", data_lines)
  ]
  
  if (length(data_lines) == 0) {
    stop("No VEP annotation records were found.")
  }
  
  text_input <- paste(
    c(header, data_lines),
    collapse = "\n"
  )
  
  read.delim(
    text = text_input,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    check.names = FALSE,
    na.strings = c("", "-", "NA")
  )
}

vep <- read_vep_table(vep_file)

required_vep_columns <- c(
  "Uploaded_variation",
  "Location",
  "Allele",
  "Gene",
  "Feature",
  "Feature_type",
  "Consequence",
  "IMPACT",
  "SYMBOL",
  "BIOTYPE",
  "CANONICAL"
)

missing_vep_columns <- setdiff(
  required_vep_columns,
  colnames(vep)
)

if (length(missing_vep_columns) > 0) {
  stop(
    "Missing VEP columns: ",
    paste(missing_vep_columns, collapse = ", ")
  )
}

# -------------------------------------------------------------------------
# Normalise fields
# -------------------------------------------------------------------------

vep$SNP <- trimws(
  as.character(vep$Uploaded_variation)
)

vep$SYMBOL <- trimws(
  as.character(vep$SYMBOL)
)

vep$Gene <- trimws(
  as.character(vep$Gene)
)

vep$Feature <- trimws(
  as.character(vep$Feature)
)

vep$Consequence <- trimws(
  as.character(vep$Consequence)
)

vep$IMPACT <- toupper(
  trimws(as.character(vep$IMPACT))
)

vep$BIOTYPE <- trimws(
  as.character(vep$BIOTYPE)
)

vep$CANONICAL <- toupper(
  trimws(as.character(vep$CANONICAL))
)

# --nearest symbol kullanıldıysa NEAREST sütunu bulunabilir.
if (!"NEAREST" %in% colnames(vep)) {
  vep$NEAREST <- NA_character_
}

vep$NEAREST <- trimws(
  as.character(vep$NEAREST)
)

# Doğrudan örtüşen gen varsa SYMBOL; yoksa nearest gene.
vep$AssignedGene <- ifelse(
  !is.na(vep$SYMBOL) & nzchar(vep$SYMBOL),
  vep$SYMBOL,
  vep$NEAREST
)

vep$GeneAssignmentType <- ifelse(
  !is.na(vep$SYMBOL) & nzchar(vep$SYMBOL),
  "overlapping_gene",
  ifelse(
    !is.na(vep$NEAREST) & nzchar(vep$NEAREST),
    "nearest_gene",
    "unassigned"
  )
)

# -------------------------------------------------------------------------
# Ranking rules
# -------------------------------------------------------------------------

impact_rank <- c(
  "HIGH" = 1L,
  "MODERATE" = 2L,
  "LOW" = 3L,
  "MODIFIER" = 4L
)

consequence_rank <- c(
  "transcript_ablation" = 1L,
  "splice_acceptor_variant" = 2L,
  "splice_donor_variant" = 3L,
  "stop_gained" = 4L,
  "frameshift_variant" = 5L,
  "stop_lost" = 6L,
  "start_lost" = 7L,
  "transcript_amplification" = 8L,
  "inframe_insertion" = 9L,
  "inframe_deletion" = 10L,
  "missense_variant" = 11L,
  "protein_altering_variant" = 12L,
  "splice_region_variant" = 13L,
  "splice_donor_region_variant" = 14L,
  "splice_polypyrimidine_tract_variant" = 15L,
  "incomplete_terminal_codon_variant" = 16L,
  "start_retained_variant" = 17L,
  "stop_retained_variant" = 18L,
  "synonymous_variant" = 19L,
  "coding_sequence_variant" = 20L,
  "mature_miRNA_variant" = 21L,
  "5_prime_UTR_variant" = 22L,
  "3_prime_UTR_variant" = 23L,
  "non_coding_transcript_exon_variant" = 24L,
  "intron_variant" = 25L,
  "NMD_transcript_variant" = 26L,
  "non_coding_transcript_variant" = 27L,
  "upstream_gene_variant" = 28L,
  "downstream_gene_variant" = 29L,
  "intergenic_variant" = 30L
)

get_consequence_rank <- function(x) {
  
  consequences <- strsplit(
    ifelse(is.na(x), "", x),
    ",",
    fixed = TRUE
  )[[1]]
  
  consequences <- trimws(consequences)
  
  ranks <- unname(
    consequence_rank[consequences]
  )
  
  ranks <- ranks[
    !is.na(ranks)
  ]
  
  if (length(ranks) == 0) {
    return(999L)
  }
  
  min(ranks)
}

vep$ImpactRank <- unname(
  impact_rank[vep$IMPACT]
)

vep$ImpactRank[
  is.na(vep$ImpactRank)
] <- 999L

vep$ConsequenceRank <- vapply(
  vep$Consequence,
  get_consequence_rank,
  integer(1)
)

vep$CanonicalRank <- ifelse(
  vep$CANONICAL == "YES",
  1L,
  2L
)

vep$ProteinCodingRank <- ifelse(
  vep$BIOTYPE == "protein_coding",
  1L,
  2L
)

vep$GeneAssignmentRank <- ifelse(
  vep$GeneAssignmentType == "overlapping_gene",
  1L,
  ifelse(
    vep$GeneAssignmentType == "nearest_gene",
    2L,
    3L
  )
)

# -------------------------------------------------------------------------
# Long SNP–gene/transcript association table
# -------------------------------------------------------------------------

long_columns <- c(
  "SNP",
  "Location",
  "Allele",
  "Gene",
  "AssignedGene",
  "GeneAssignmentType",
  "Feature",
  "Feature_type",
  "Consequence",
  "IMPACT",
  "BIOTYPE",
  "CANONICAL",
  "VARIANT_CLASS",
  "DISTANCE",
  "Amino_acids",
  "Codons",
  "Protein_position",
  "HGNC_ID",
  "NEAREST"
)

long_columns <- intersect(
  long_columns,
  colnames(vep)
)

long_table <- unique(
  vep[, long_columns, drop = FALSE]
)

write.table(
  long_table,
  file = long_output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  na = ""
)

# -------------------------------------------------------------------------
# Select one representative annotation per SNP
#
# Priority:
#   1. Overlapping gene rather than nearest gene
#   2. More severe IMPACT
#   3. More severe consequence
#   4. Canonical transcript
#   5. Protein-coding transcript
# -------------------------------------------------------------------------

vep_ordered <- vep[
  order(
    vep$SNP,
    vep$GeneAssignmentRank,
    vep$ImpactRank,
    vep$ConsequenceRank,
    vep$CanonicalRank,
    vep$ProteinCodingRank,
    vep$AssignedGene,
    vep$Feature,
    na.last = TRUE
  ),
  ,
  drop = FALSE
]

primary <- vep_ordered[
  !duplicated(vep_ordered$SNP),
  ,
  drop = FALSE
]

# Her SNP için tüm örtüşen genleri de koru.
gene_list <- aggregate(
  AssignedGene ~ SNP,
  data = vep[
    !is.na(vep$AssignedGene) &
      nzchar(vep$AssignedGene),
    c("SNP", "AssignedGene")
  ],
  FUN = function(x) {
    paste(
      sort(unique(x)),
      collapse = ";"
    )
  }
)

colnames(gene_list)[2] <- "AllAssignedGenes"

transcript_count <- aggregate(
  Feature ~ SNP,
  data = vep[
    !is.na(vep$Feature) &
      nzchar(vep$Feature),
    c("SNP", "Feature")
  ],
  FUN = function(x) {
    length(unique(x))
  }
)

colnames(transcript_count)[2] <- "TranscriptCount"

gene_count <- aggregate(
  AssignedGene ~ SNP,
  data = vep[
    !is.na(vep$AssignedGene) &
      nzchar(vep$AssignedGene),
    c("SNP", "AssignedGene")
  ],
  FUN = function(x) {
    length(unique(x))
  }
)

colnames(gene_count)[2] <- "AssignedGeneCount"

summary_columns <- c(
  "SNP",
  "Location",
  "Allele",
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
  "DISTANCE",
  "Amino_acids",
  "Codons",
  "Protein_position",
  "HGNC_ID",
  "NEAREST"
)

summary_columns <- intersect(
  summary_columns,
  colnames(primary)
)

summary_table <- primary[
  ,
  summary_columns,
  drop = FALSE
]

summary_table <- merge(
  summary_table,
  gene_list,
  by = "SNP",
  all.x = TRUE,
  sort = FALSE
)

summary_table <- merge(
  summary_table,
  gene_count,
  by = "SNP",
  all.x = TRUE,
  sort = FALSE
)

summary_table <- merge(
  summary_table,
  transcript_count,
  by = "SNP",
  all.x = TRUE,
  sort = FALSE
)

summary_table <- summary_table[
  match(
    unique(primary$SNP),
    summary_table$SNP
  ),
  ,
  drop = FALSE
]

write.table(
  summary_table,
  file = summary_output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  na = ""
)

# -------------------------------------------------------------------------
# Merge with InstrumentEvidence
# -------------------------------------------------------------------------

instrument <- read.delim(
  instrument_file,
  sep = "\t",
  header = TRUE,
  quote = "",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (!"SNP" %in% colnames(instrument)) {
  stop("InstrumentEvidence.tsv does not contain an SNP column.")
}

instrument$SNP <- trimws(
  as.character(instrument$SNP)
)

# Her orijinal instrument satırına benzersiz sıra numarası ver.
instrument$InstrumentRowID <- seq_len(
  nrow(instrument)
)

# summary_table her SNP için yalnızca bir satır içermeli.
if (anyDuplicated(summary_table$SNP) > 0) {
  duplicated_snps <- unique(
    summary_table$SNP[
      duplicated(summary_table$SNP) |
        duplicated(summary_table$SNP, fromLast = TRUE)
    ]
  )
  
  stop(
    "VariantGeneSummary contains duplicated SNPs: ",
    paste(head(duplicated_snps, 20), collapse = ", ")
  )
}

merged <- merge(
  instrument,
  summary_table,
  by = "SNP",
  all.x = TRUE,
  sort = FALSE
)

# Orijinal InstrumentEvidence sırasını geri getir.
merged <- merged[
  order(merged$InstrumentRowID),
  ,
  drop = FALSE
]

merged$InstrumentRowID <- NULL

write.table(
  merged,
  file = merged_output,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  na = ""
)
# -------------------------------------------------------------------------
# Report
# -------------------------------------------------------------------------

cat("\nVEP annotation parsing completed.\n\n")

cat("Input files:\n")
cat("  VEP:", vep_file, "\n")
cat("  Instrument evidence:", instrument_file, "\n\n")

cat("VEP transcript records:", nrow(vep), "\n")
cat("Unique annotated SNPs:", length(unique(vep$SNP)), "\n")
cat(
  "Unique assigned genes:",
  length(unique(
    vep$AssignedGene[
      !is.na(vep$AssignedGene) &
        nzchar(vep$AssignedGene)
    ]
  )),
  "\n"
)

cat(
  "Missense primary annotations:",
  sum(
    grepl(
      "missense_variant",
      summary_table$Consequence,
      fixed = TRUE
    ),
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Instrument rows with assigned gene:",
  sum(
    !is.na(merged$AssignedGene) &
      nzchar(merged$AssignedGene)
  ),
  "of",
  nrow(merged),
  "\n\n"
)

cat("Outputs:\n")
cat("  ", long_output, "\n")
cat("  ", summary_output, "\n")
cat("  ", merged_output, "\n")
