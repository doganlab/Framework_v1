#!/usr/bin/env Rscript

# ============================================================
# Framework_v1
# 18_post_prioritization.R
#
# Purpose
# -------
# Build Supplementary_Data_S6.xlsx as a POST-PRIORITIZATION,
# CAD-SPECIFIC functional / V2G support layer by integrating:
#
#   Framework Supplementary_Data_S5.xlsx
#     - AllModules
#     - TierA / TierB / TierC
#
# with Aragam et al. 2022 CAD supplementary tables:
#     ST20  functionally informed fine-mapping regions
#     ST21  95% credible-set variants
#     ST22  high-PPA variants + GTEx/STARNET eQTL support
#     ST27  GTEx eQTLs in high LD with CAD sentinel variants
#     ST28  STARNET eQTLs in high LD with CAD sentinel variants
#     ST31  integrated eight-predictor CAD causal-gene prioritization
#
# Design principles
# -----------------
# 1) This script DOES NOT change the Framework_v1 candidate tiers.
# 2) Aragam evidence is added only after framework prioritization.
# 3) Direct variant/locus evidence is kept separate from gene-level
#    CAD support.
# 4) ST27/ST28 are described as high-LD eQTL support, NOT formal
#    eQTL colocalization.
# 5) Absence from an Aragam supplementary table is NOT interpreted
#    as evidence against a framework candidate.
# 6) No composite biological score is calculated. "N_AragamSupportAxes"
#    is a descriptive count of available evidence axes only.
#
# Expected location when saved in the project:
#   Framework_v1/scripts/29_build_Aragam_S6.R
#
# Output:
#   results/candidate_prioritization/Aragam_CAD_support/
#     Supplementary_Data_S6.xlsx
# ============================================================

options(stringsAsFactors = FALSE, warn = 1, scipen = 999)

# ------------------------------------------------------------
# Bootstrap / project paths
# ------------------------------------------------------------

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
if (is.na(script_arg) || !nzchar(script_arg)) {
  stop("This script must be run with Rscript.", call. = FALSE)
}

script_path <- normalizePath(
  sub("^--file=", "", script_arg),
  winslash = "/",
  mustWork = TRUE
)

source(file.path(dirname(script_path), "..", "R", "bootstrap.R"))

required_packages <- c("data.table", "readxl", "openxlsx")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Missing R package(s): ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them before running this script.",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(openxlsx)
})

ARAGAM_DIR <- file.path(
  DIR_RAW,
  "AragamSupplementaryTables"
)

ARAGAM_FILE <- file.path(
  ARAGAM_DIR,
  "41588_2022_1233_MOESM4_ESM.xlsx"
)

S5_FILE <- file.path(
  ARAGAM_DIR,
  "Supplementary_Data_S5.xlsx"
)

OUTPUT_DIR <- file.path(
  DIR_RESULTS,
  "candidate_prioritization",
  "Aragam_CAD_support"
)

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUTPUT_FILE <- file.path(
  OUTPUT_DIR,
  "Supplementary_Data_S6.xlsx"
)

QC_TSV <- file.path(
  OUTPUT_DIR,
  "Supplementary_Data_S6_QC.tsv"
)

for (f in c(ARAGAM_FILE, S5_FILE)) {
  if (!file.exists(f)) {
    stop("Required input file not found: ", f, call. = FALSE)
  }
}

cat("\n=== Framework_v1: Aragam CAD-specific support layer ===\n")
cat("Aragam workbook: ", ARAGAM_FILE, "\n", sep = "")
cat("Framework S5:    ", S5_FILE, "\n", sep = "")
cat("Output:          ", OUTPUT_FILE, "\n\n", sep = "")

# ============================================================
# Helpers
# ============================================================

clean_name <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("[[:space:]]+", "_", x)
  x <- gsub("[^A-Za-z0-9_.]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

norm_gene <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x %in% c("", "NA", "N/A", "-")] <- NA_character_
  x
}

norm_snp <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x[!grepl("^rs[0-9]+$", x)] <- NA_character_
  x
}

collapse_unique <- function(x, sep = ";") {
  x <- as.character(x)
  x <- trimws(x)
  x <- sort(unique(x[!is.na(x) & x != "" & x != "NA"]))
  if (!length(x)) return(NA_character_)
  paste(x, collapse = sep)
}

safe_max_num <- function(x) {
  z <- suppressWarnings(as.numeric(x))
  if (all(is.na(z))) return(NA_real_)
  max(z, na.rm = TRUE)
}

safe_min_num <- function(x) {
  z <- suppressWarnings(as.numeric(x))
  if (all(is.na(z))) return(NA_real_)
  min(z, na.rm = TRUE)
}

safe_max_int <- function(x) {
  z <- suppressWarnings(as.integer(x))
  if (all(is.na(z))) return(NA_integer_)
  max(z, na.rm = TRUE)
}

safe_min_int <- function(x) {
  z <- suppressWarnings(as.integer(x))
  if (all(is.na(z))) return(NA_integer_)
  min(z, na.rm = TRUE)
}

as_flag <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  toupper(trimws(as.character(x))) %in% c("TRUE", "T", "YES", "Y", "1")
}

# Split strings that can contain one or more gene symbols.
# Used only for exact gene-token matching, never fuzzy matching.
gene_tokens <- function(x) {
  if (length(x) == 0L || all(is.na(x))) return(character())
  y <- paste(as.character(x[!is.na(x)]), collapse = ";")
  z <- unlist(strsplit(y, "[^A-Za-z0-9.-]+", perl = TRUE), use.names = FALSE)
  z <- norm_gene(z)
  sort(unique(z[!is.na(z) & z != ""]))
}

contains_gene <- function(x, gene) {
  gene <- norm_gene(gene)
  if (is.na(gene) || !nzchar(gene)) return(FALSE)
  gene %in% gene_tokens(x)
}

parse_rsids <- function(x) {
  x <- paste(as.character(x[!is.na(x)]), collapse = ";")
  hits <- regmatches(
    x,
    gregexpr("rs[0-9]+", x, ignore.case = TRUE, perl = TRUE)
  )[[1]]
  if (!length(hits) || identical(hits, character(0))) return(character())
  sort(unique(tolower(hits)))
}

resolve_col <- function(dt, aliases, required = FALSE, label = NULL) {
  hit <- aliases[aliases %in% names(dt)]
  if (length(hit)) return(hit[1])

  # Also try normalized names because Excel punctuation can vary.
  norm_names <- clean_name(names(dt))
  norm_alias <- clean_name(aliases)
  idx <- match(norm_alias, norm_names, nomatch = 0L)
  idx <- idx[idx > 0L]
  if (length(idx)) return(names(dt)[idx[1]])

  if (required) {
    if (is.null(label)) label <- paste(aliases, collapse = " / ")
    stop(
      "Required column not found: ", label,
      "\nAvailable columns: ", paste(names(dt), collapse = ", "),
      call. = FALSE
    )
  }
  NA_character_
}

get_col <- function(dt, aliases, default = NA_character_) {
  nm <- resolve_col(dt, aliases, required = FALSE)
  if (is.na(nm)) return(rep(default, nrow(dt)))
  dt[[nm]]
}

# Read an Excel sheet whose actual column header may start after title/note rows.
# The first row containing at least one expected header token is used as header.
read_sheet_detect_header <- function(path, sheet, expected_tokens) {
  raw <- readxl::read_excel(
    path,
    sheet = sheet,
    col_names = FALSE,
    col_types = "text",
    .name_repair = "minimal"
  )

  if (!nrow(raw)) stop("Sheet is empty: ", sheet, call. = FALSE)

  token_norm <- clean_name(expected_tokens)
  hit_row <- NA_integer_

  for (i in seq_len(nrow(raw))) {
    vals <- clean_name(unlist(raw[i, ], use.names = FALSE))
    if (any(vals %in% token_norm, na.rm = TRUE)) {
      hit_row <- i
      break
    }
  }

  if (is.na(hit_row)) {
    stop(
      "Could not detect header row in ", sheet,
      ". Expected one of: ", paste(expected_tokens, collapse = ", "),
      call. = FALSE
    )
  }

  x <- data.table::as.data.table(
    readxl::read_excel(
      path,
      sheet = sheet,
      skip = hit_row - 1L,
      col_names = TRUE,
      .name_repair = "unique"
    )
  )

  # Drop completely empty rows/columns.
  if (nrow(x)) {
    keep_row <- apply(
      as.data.frame(x),
      1,
      function(v) any(!is.na(v) & trimws(as.character(v)) != "")
    )
    x <- x[keep_row]
  }

  empty_cols <- vapply(
    x,
    function(v) all(is.na(v) | trimws(as.character(v)) == ""),
    logical(1)
  )
  if (any(empty_cols)) x[, names(x)[empty_cols] := NULL]

  x[]
}

# ============================================================
# 1. Load Framework S5
# ============================================================

s5_sheets <- readxl::excel_sheets(S5_FILE)
required_s5_sheets <- c("AllModules", "TierA", "TierB", "TierC")
missing_s5 <- setdiff(required_s5_sheets, s5_sheets)

if (length(missing_s5)) {
  stop(
    "Supplementary_Data_S5.xlsx is missing sheet(s): ",
    paste(missing_s5, collapse = ", "),
    call. = FALSE
  )
}

read_s5_sheet <- function(sheet) {
  data.table::as.data.table(
    readxl::read_excel(
      S5_FILE,
      sheet = sheet,
      .name_repair = "unique"
    )
  )
}

all_modules <- read_s5_sheet("AllModules")
tierA_s5 <- read_s5_sheet("TierA")
tierB_s5 <- read_s5_sheet("TierB")
tierC_s5 <- read_s5_sheet("TierC")

required_module_cols <- c("Gene", "NatMed_DrugCategory", "SNPs")
miss <- setdiff(required_module_cols, names(all_modules))
if (length(miss)) {
  stop(
    "S5 AllModules is missing required columns: ",
    paste(miss, collapse = ", "),
    call. = FALSE
  )
}

all_modules[, Gene := norm_gene(Gene)]
all_modules[, S5_Row := .I]
all_modules[, FrameworkModuleID := paste0(
  "M", sprintf("%05d", S5_Row), "::",
  Gene, "::", NatMed_DrugCategory
)]

# Preserve an existing index if available; otherwise S5_Row is the stable identifier.
index_col <- resolve_col(all_modules, c("index", "Index", "S5_Index"), required = FALSE)
if (!is.na(index_col)) {
  all_modules[, S5_Index := as.character(get(index_col))]
} else {
  all_modules[, S5_Index := as.character(S5_Row)]
}

# Expand unique framework SNPs once. This avoids duplicating Aragam lookup work
# for the same rsID across many gene-drug modules.
all_framework_snps <- sort(unique(unlist(lapply(all_modules$SNPs, parse_rsids))))
all_framework_snps <- all_framework_snps[!is.na(all_framework_snps)]

framework_genes <- sort(unique(all_modules$Gene[!is.na(all_modules$Gene)]))

cat("S5 modules:       ", nrow(all_modules), "\n", sep = "")
cat("Unique S5 genes:  ", length(framework_genes), "\n", sep = "")
cat("Unique S5 rsIDs:  ", length(all_framework_snps), "\n\n", sep = "")

# ============================================================
# 2. Load Aragam sheets
# ============================================================

arag_sheets <- readxl::excel_sheets(ARAGAM_FILE)
needed_arag <- c("ST20", "ST21", "ST22", "ST27", "ST28", "ST31")
missing_arag <- setdiff(needed_arag, arag_sheets)

if (length(missing_arag)) {
  stop(
    "Aragam workbook is missing sheet(s): ",
    paste(missing_arag, collapse = ", "),
    call. = FALSE
  )
}

st20 <- read_sheet_detect_header(
  ARAGAM_FILE, "ST20",
  c("FGWAS_region", "Lead_variant_rsID", "BestSNP_PPA")
)

st21 <- read_sheet_detect_header(
  ARAGAM_FILE, "ST21",
  c("FGWAS_region", "rsID", "PPA")
)

st22 <- read_sheet_detect_header(
  ARAGAM_FILE, "ST22",
  c("FGWAS_region", "FGWAS_lead_variant", "FGWAS_proxy_SNP", "PPA")
)

st27 <- read_sheet_detect_header(
  ARAGAM_FILE, "ST27",
  c("Sentinel_rsID", "Proxy_rsID", "Gene_name", "Tissue")
)

st28 <- read_sheet_detect_header(
  ARAGAM_FILE, "ST28",
  c("CAD_variant_rsID", "Index_eSNP", "Gene_symbol", "Tissue")
)

st31 <- read_sheet_detect_header(
  ARAGAM_FILE, "ST31",
  c("rsID", "most_likely_causal_gene", "n_predictors", "95%_credible_set")
)

# ============================================================
# 3. Standardize Aragam lookup tables
# ============================================================

# ---------- ST20: region-level fine-mapping ----------
st20_std <- data.table(
  FGWAS_region = as.character(get_col(st20, c("FGWAS_region"))),
  ST20_LeadVariant = norm_snp(get_col(st20, c("Lead_variant_rsID", "Lead_variant", "lead_variant_rsID"))),
  ST20_LocusName = as.character(get_col(st20, c("LocusName", "Locus_Name", "locus_name"))),
  ST20_BestSNP_PPA = suppressWarnings(as.numeric(get_col(st20, c("BestSNP_PPA"), default = NA_real_))),
  ST20_CredibleSetSize = suppressWarnings(as.integer(get_col(
    st20,
    c("Variants_in_95%_credible_set", "95%_credible_set", "Variants_in_95_credible_set"),
    default = NA_integer_
  ))),
  ST20_StrongestTissue = as.character(get_col(
    st20,
    c("Most_enriched_tissue", "strongest_tissue", "Strongest_tissue")
  ))
)

st20_std <- st20_std[!is.na(FGWAS_region) & FGWAS_region != ""]
st20_std <- unique(st20_std, by = "FGWAS_region")

# ---------- ST21: 95% credible-set variants ----------
st21_std <- data.table(
  SNP = norm_snp(get_col(st21, c("rsID", "rsid"))),
  FGWAS_region = as.character(get_col(st21, c("FGWAS_region"))),
  ST21_PPA = suppressWarnings(as.numeric(get_col(st21, c("PPA"), default = NA_real_)))
)
st21_std <- st21_std[!is.na(SNP)]

st21_std <- merge(
  st21_std,
  st20_std,
  by = "FGWAS_region",
  all.x = TRUE,
  sort = FALSE
)

st21_snp <- st21_std[, .(
  ST21_CredibleSet = TRUE,
  ST21_MaxPPA = safe_max_num(ST21_PPA),
  ST21_Regions = collapse_unique(FGWAS_region),
  ST21_Loci = collapse_unique(ST20_LocusName),
  ST20_MaxBestSNP_PPA = safe_max_num(ST20_BestSNP_PPA),
  ST20_MinCredibleSetSize = safe_min_int(ST20_CredibleSetSize),
  ST20_StrongestTissues = collapse_unique(ST20_StrongestTissue)
), by = SNP]

# ---------- ST22: high-PPA likely causal variants + eQTL context ----------
st22_lead_col <- resolve_col(
  st22,
  c("FGWAS_lead_variant", "FGWAS_lead_variant_rsID", "Lead_variant_rsID"),
  required = FALSE
)
st22_proxy_col <- resolve_col(
  st22,
  c("FGWAS_proxy_SNP", "Proxy_SNP", "proxy_SNP"),
  required = FALSE
)

st22_base <- copy(st22)
st22_base[, ST22_RowID := .I]
st22_base[, ST22_PPA := suppressWarnings(as.numeric(get_col(st22_base, c("PPA"), default = NA_real_)))]
st22_base[, ST22_Region := as.character(get_col(st22_base, c("FGWAS_region")))]
st22_base[, ST22_Locus := as.character(get_col(st22_base, c("locus_name", "LocusName", "locus")))]

st22_base[, ST22_GTExGene := norm_gene(get_col(
  st22_base,
  c("gtex_gene", "GTEx_gene", "GTEx_gene_name", "gtex_gene_name")
))]
st22_base[, ST22_GTExTissue := as.character(get_col(
  st22_base,
  c("gtex_tissue", "GTEx_tissue")
))]
st22_base[, ST22_GTExR2 := suppressWarnings(as.numeric(get_col(
  st22_base,
  c("gtex_r2", "GTEx_r2"), default = NA_real_
)))]

st22_base[, ST22_STARNETGene := norm_gene(get_col(
  st22_base,
  c("starnet_gene", "STARNET_gene", "STARNET_gene_name", "starnet_gene_name")
))]
st22_base[, ST22_STARNETTissue := as.character(get_col(
  st22_base,
  c("starnet_tissue", "STARNET_tissue")
))]
st22_base[, ST22_STARNETR2 := suppressWarnings(as.numeric(get_col(
  st22_base,
  c("starnet_r2", "STARNET_r2"), default = NA_real_
)))]

st22_long_parts <- list()
if (!is.na(st22_lead_col)) {
  st22_long_parts[["lead"]] <- st22_base[, .(
    SNP = norm_snp(get(st22_lead_col)),
    ST22_MatchRole = "FGWAS_lead_variant",
    ST22_RowID,
    ST22_PPA,
    ST22_Region,
    ST22_Locus,
    ST22_GTExGene,
    ST22_GTExTissue,
    ST22_GTExR2,
    ST22_STARNETGene,
    ST22_STARNETTissue,
    ST22_STARNETR2
  )]
}
if (!is.na(st22_proxy_col)) {
  st22_long_parts[["proxy"]] <- st22_base[, .(
    SNP = norm_snp(get(st22_proxy_col)),
    ST22_MatchRole = "FGWAS_proxy_SNP",
    ST22_RowID,
    ST22_PPA,
    ST22_Region,
    ST22_Locus,
    ST22_GTExGene,
    ST22_GTExTissue,
    ST22_GTExR2,
    ST22_STARNETGene,
    ST22_STARNETTissue,
    ST22_STARNETR2
  )]
}

st22_long <- rbindlist(st22_long_parts, fill = TRUE)
st22_long <- st22_long[!is.na(SNP)]

st22_snp <- st22_long[, .(
  ST22_HighPPA_TableMatch = TRUE,
  ST22_MaxPPA = safe_max_num(ST22_PPA),
  ST22_MatchRoles = collapse_unique(ST22_MatchRole),
  ST22_Regions = collapse_unique(ST22_Region),
  ST22_Loci = collapse_unique(ST22_Locus),
  ST22_GTExGenes = collapse_unique(ST22_GTExGene),
  ST22_GTExTissues = collapse_unique(ST22_GTExTissue),
  ST22_MaxGTExR2 = safe_max_num(ST22_GTExR2),
  ST22_STARNETGenes = collapse_unique(ST22_STARNETGene),
  ST22_STARNETTissues = collapse_unique(ST22_STARNETTissue),
  ST22_MaxSTARNETR2 = safe_max_num(ST22_STARNETR2)
), by = SNP]

# ---------- ST27: GTEx eQTLs in high LD with CAD sentinels ----------
st27_std <- data.table(
  ST27_SentinelSNP = norm_snp(get_col(st27, c("Sentinel_rsID"))),
  ST27_ProxySNP = norm_snp(get_col(st27, c("Proxy_rsID"))),
  ST27_r2 = suppressWarnings(as.numeric(get_col(st27, c("r2", "R2"), default = NA_real_))),
  ST27_Gene = norm_gene(get_col(st27, c("Gene_name", "Gene", "gene_name"))),
  ST27_Tissue = as.character(get_col(st27, c("Tissue", "tissue"))),
  ST27_Slope = suppressWarnings(as.numeric(get_col(st27, c("Slope", "beta", "Beta"), default = NA_real_))),
  ST27_P = suppressWarnings(as.numeric(get_col(st27, c("Nominal_p-value", "Nominal_p_value", "P-value", "P_value"), default = NA_real_)))
)

st27_long <- rbindlist(list(
  st27_std[!is.na(ST27_SentinelSNP), .(
    SNP = ST27_SentinelSNP,
    ST27_MatchRole = "CAD_sentinel",
    ST27_SentinelSNP,
    ST27_ProxySNP,
    ST27_r2,
    ST27_Gene,
    ST27_Tissue,
    ST27_Slope,
    ST27_P
  )],
  st27_std[!is.na(ST27_ProxySNP), .(
    SNP = ST27_ProxySNP,
    ST27_MatchRole = "GTEx_eQTL_proxy",
    ST27_SentinelSNP,
    ST27_ProxySNP,
    ST27_r2,
    ST27_Gene,
    ST27_Tissue,
    ST27_Slope,
    ST27_P
  )]
), fill = TRUE)

st27_snp <- st27_long[, .(
  ST27_GTEx_HighLD_Support = TRUE,
  ST27_MatchRoles = collapse_unique(ST27_MatchRole),
  ST27_CADSentinels = collapse_unique(ST27_SentinelSNP),
  ST27_ProxySNPs = collapse_unique(ST27_ProxySNP),
  ST27_MaxR2 = safe_max_num(ST27_r2),
  ST27_GTExGenes = collapse_unique(ST27_Gene),
  ST27_GTExTissues = collapse_unique(ST27_Tissue),
  ST27_MinP = safe_min_num(ST27_P)
), by = SNP]

# ---------- ST28: STARNET eQTLs in high LD with CAD sentinels ----------
st28_std <- data.table(
  ST28_CADVariant = norm_snp(get_col(st28, c("CAD_variant_rsID", "CAD_variant"))),
  ST28_IndexeSNP = norm_snp(get_col(st28, c("Index_eSNP", "Index_eSNP_rsID"))),
  ST28_r2 = suppressWarnings(as.numeric(get_col(st28, c("r2", "R2"), default = NA_real_))),
  ST28_Gene = norm_gene(get_col(st28, c("Gene_symbol", "Gene", "gene_symbol"))),
  ST28_Tissue = as.character(get_col(st28, c("Tissue", "tissue"))),
  ST28_Beta = suppressWarnings(as.numeric(get_col(st28, c("Beta", "beta"), default = NA_real_))),
  ST28_P = suppressWarnings(as.numeric(get_col(st28, c("P-value", "P_value", "Pvalue"), default = NA_real_)))
)

st28_long <- rbindlist(list(
  st28_std[!is.na(ST28_CADVariant), .(
    SNP = ST28_CADVariant,
    ST28_MatchRole = "CAD_sentinel",
    ST28_CADVariant,
    ST28_IndexeSNP,
    ST28_r2,
    ST28_Gene,
    ST28_Tissue,
    ST28_Beta,
    ST28_P
  )],
  st28_std[!is.na(ST28_IndexeSNP), .(
    SNP = ST28_IndexeSNP,
    ST28_MatchRole = "STARNET_index_eSNP",
    ST28_CADVariant,
    ST28_IndexeSNP,
    ST28_r2,
    ST28_Gene,
    ST28_Tissue,
    ST28_Beta,
    ST28_P
  )]
), fill = TRUE)

st28_snp <- st28_long[, .(
  ST28_STARNET_HighLD_Support = TRUE,
  ST28_MatchRoles = collapse_unique(ST28_MatchRole),
  ST28_CADSentinels = collapse_unique(ST28_CADVariant),
  ST28_IndexeSNPs = collapse_unique(ST28_IndexeSNP),
  ST28_MaxR2 = safe_max_num(ST28_r2),
  ST28_STARNETGenes = collapse_unique(ST28_Gene),
  ST28_STARNETTissues = collapse_unique(ST28_Tissue),
  ST28_MinP = safe_min_num(ST28_P)
), by = SNP]

# ---------- ST31: integrated causal-gene prioritization ----------
st31_std <- data.table(
  SNP = norm_snp(get_col(st31, c("rsID", "rsid"))),
  ST31_Locus = as.character(get_col(st31, c("Locus", "locus"))),
  ST31_Marker = as.character(get_col(st31, c("markername", "MarkerName", "marker_name"))),
  ST31_CredibleSetSize = suppressWarnings(as.integer(get_col(
    st31,
    c("95%_credible_set", "95_credible_set"),
    default = NA_integer_
  ))),
  ST31_BestSNP_PPA = suppressWarnings(as.numeric(get_col(st31, c("BestSNP_PPA"), default = NA_real_))),
  ST31_NearestGene = norm_gene(get_col(st31, c("nearest_gene"))),
  ST31_RareVariantGene = norm_gene(get_col(st31, c("previous_rare_variant"))),
  ST31_MonogenicGene = norm_gene(get_col(st31, c("monogenic_disorder"))),
  ST31_DrugMRGene = norm_gene(get_col(st31, c("drug/MR_evidence", "drug_MR_evidence", "drug.MR_evidence"))),
  ST31_ProteinAlteringGene = norm_gene(get_col(st31, c("protein_altering_variant"))),
  ST31_PrimaryPoPSGene = norm_gene(get_col(st31, c("Primary_PoPS_gene"))),
  ST31_SecondaryPoPSGene = norm_gene(get_col(st31, c("Secondary_PoPS_gene"))),
  ST31_TopGTExGene = norm_gene(get_col(st31, c("top_GTEx_eQTL_gene"))),
  ST31_OtherGTExGenes = as.character(get_col(st31, c("other_GTEx_eQTL_genes"))),
  ST31_TopSTARNETGene = norm_gene(get_col(st31, c("top_STARNET_eQTL_gene"))),
  ST31_OtherSTARNETGenes = as.character(get_col(st31, c("other_STARNET_eQTL_genes"))),
  ST31_MouseGene = norm_gene(get_col(st31, c("mouse_phenotype"))),
  ST31_MostLikelyCausalGene = norm_gene(get_col(st31, c("most_likely_causal_gene"))),
  ST31_NPredictors = suppressWarnings(as.integer(get_col(st31, c("n_predictors"), default = NA_integer_)))
)

st31_std <- st31_std[!is.na(SNP)]

st31_snp <- st31_std[, .(
  ST31_CADSentinelMatch = TRUE,
  ST31_Loci = collapse_unique(ST31_Locus),
  ST31_Markers = collapse_unique(ST31_Marker),
  ST31_MinCredibleSetSize = safe_min_int(ST31_CredibleSetSize),
  ST31_MaxBestSNP_PPA = safe_max_num(ST31_BestSNP_PPA),
  ST31_MostLikelyCausalGenes = collapse_unique(ST31_MostLikelyCausalGene),
  ST31_MaxNPredictors = safe_max_int(ST31_NPredictors),
  ST31_NearestGenes = collapse_unique(ST31_NearestGene),
  ST31_PrimaryPoPSGenes = collapse_unique(ST31_PrimaryPoPSGene),
  ST31_SecondaryPoPSGenes = collapse_unique(ST31_SecondaryPoPSGene),
  ST31_TopGTExGenes = collapse_unique(ST31_TopGTExGene),
  ST31_OtherGTExGenes = collapse_unique(ST31_OtherGTExGenes),
  ST31_TopSTARNETGenes = collapse_unique(ST31_TopSTARNETGene),
  ST31_OtherSTARNETGenes = collapse_unique(ST31_OtherSTARNETGenes),
  ST31_ProteinAlteringGenes = collapse_unique(ST31_ProteinAlteringGene),
  ST31_DrugMRGenes = collapse_unique(ST31_DrugMRGene),
  ST31_RareVariantGenes = collapse_unique(ST31_RareVariantGene),
  ST31_MonogenicGenes = collapse_unique(ST31_MonogenicGene),
  ST31_MouseGenes = collapse_unique(ST31_MouseGene)
), by = SNP]

# ============================================================
# 4. Build unique-SNP Aragam support table
# ============================================================

variant_support <- data.table(SNP = all_framework_snps)

for (lookup in list(st21_snp, st22_snp, st27_snp, st28_snp, st31_snp)) {
  variant_support <- merge(
    variant_support,
    lookup,
    by = "SNP",
    all.x = TRUE,
    sort = FALSE
  )
}

flag_names <- c(
  "ST21_CredibleSet",
  "ST22_HighPPA_TableMatch",
  "ST27_GTEx_HighLD_Support",
  "ST28_STARNET_HighLD_Support",
  "ST31_CADSentinelMatch"
)

for (nm in flag_names) {
  if (!nm %in% names(variant_support)) variant_support[, (nm) := FALSE]
  variant_support[is.na(get(nm)), (nm) := FALSE]
}

variant_support[, Any_FineMappingSupport :=
  ST21_CredibleSet | ST22_HighPPA_TableMatch]

variant_support[, Any_HighLD_eQTLSupport :=
  ST27_GTEx_HighLD_Support | ST28_STARNET_HighLD_Support]

variant_support[, Any_AragamVariantSupport :=
  Any_FineMappingSupport |
    Any_HighLD_eQTLSupport |
    ST31_CADSentinelMatch]

variant_support[, N_AragamSupportAxes := rowSums(cbind(
  Any_FineMappingSupport,
  ST27_GTEx_HighLD_Support,
  ST28_STARNET_HighLD_Support,
  ST31_CADSentinelMatch
), na.rm = TRUE)]

variant_support[, AragamVariantSupportClass := fifelse(
  ST21_CredibleSet | ST22_HighPPA_TableMatch,
  "Fine-mapping-supported CAD locus variant",
  fifelse(
    ST27_GTEx_HighLD_Support | ST28_STARNET_HighLD_Support,
    "CAD-sentinel/high-LD eQTL-linked variant",
    fifelse(
      ST31_CADSentinelMatch,
      "Aragam CAD sentinel variant",
      "No direct Aragam variant-table match"
    )
  )
)]

setorder(
  variant_support,
  -Any_AragamVariantSupport,
  -N_AragamSupportAxes,
  SNP
)

# ============================================================
# 5. Build gene-level Aragam support table
# ============================================================

# Gene-level support is intentionally kept separate from direct variant matching.
# It asks whether the Framework-assigned gene appears anywhere in the Aragam
# CAD-specific functional / causal-gene resources.

st27_gene_set <- sort(unique(st27_std$ST27_Gene[!is.na(st27_std$ST27_Gene)]))
st28_gene_set <- sort(unique(st28_std$ST28_Gene[!is.na(st28_std$ST28_Gene)]))
st22_gtex_gene_set <- sort(unique(st22_base$ST22_GTExGene[!is.na(st22_base$ST22_GTExGene)]))
st22_starnet_gene_set <- sort(unique(st22_base$ST22_STARNETGene[!is.na(st22_base$ST22_STARNETGene)]))

st31_mlc_set <- gene_tokens(st31_std$ST31_MostLikelyCausalGene)
st31_pops1_set <- gene_tokens(st31_std$ST31_PrimaryPoPSGene)
st31_pops2_set <- gene_tokens(st31_std$ST31_SecondaryPoPSGene)
st31_nearest_set <- gene_tokens(st31_std$ST31_NearestGene)
st31_protein_set <- gene_tokens(st31_std$ST31_ProteinAlteringGene)
st31_drugmr_set <- gene_tokens(st31_std$ST31_DrugMRGene)

# top + "other" eQTL genes need tokenization.
st31_gtex_tokens <- gene_tokens(c(st31_std$ST31_TopGTExGene, st31_std$ST31_OtherGTExGenes))
st31_starnet_tokens <- gene_tokens(c(st31_std$ST31_TopSTARNETGene, st31_std$ST31_OtherSTARNETGenes))

gene_support <- data.table(Gene = framework_genes)

gene_support[, ST31_MostLikelyCausalGene := Gene %in% st31_mlc_set]
gene_support[, ST31_PrimaryPoPS := Gene %in% st31_pops1_set]
gene_support[, ST31_SecondaryPoPS := Gene %in% st31_pops2_set]
gene_support[, ST31_NearestGene := Gene %in% st31_nearest_set]
gene_support[, ST31_ProteinAltering := Gene %in% st31_protein_set]
gene_support[, ST31_DrugMREvidence := Gene %in% st31_drugmr_set]
gene_support[, ST31_GTEx_eQTL := Gene %in% st31_gtex_tokens]
gene_support[, ST31_STARNET_eQTL := Gene %in% st31_starnet_tokens]
gene_support[, ST27_GTEx_HighLD_eQTL := Gene %in% st27_gene_set]
gene_support[, ST28_STARNET_HighLD_eQTL := Gene %in% st28_gene_set]
gene_support[, ST22_GTEx_HighPPA_eQTL := Gene %in% st22_gtex_gene_set]
gene_support[, ST22_STARNET_HighPPA_eQTL := Gene %in% st22_starnet_gene_set]

gene_support[, Any_AragamGeneSupport := rowSums(cbind(
  ST31_MostLikelyCausalGene,
  ST31_PrimaryPoPS,
  ST31_SecondaryPoPS,
  ST31_ProteinAltering,
  ST31_DrugMREvidence,
  ST31_GTEx_eQTL,
  ST31_STARNET_eQTL,
  ST27_GTEx_HighLD_eQTL,
  ST28_STARNET_HighLD_eQTL,
  ST22_GTEx_HighPPA_eQTL,
  ST22_STARNET_HighPPA_eQTL
), na.rm = TRUE) > 0]

gene_support[, N_AragamGeneEvidenceAxes := rowSums(cbind(
  ST31_MostLikelyCausalGene,
  ST31_PrimaryPoPS,
  ST31_SecondaryPoPS,
  ST31_ProteinAltering,
  ST31_DrugMREvidence,
  ST31_GTEx_eQTL | ST27_GTEx_HighLD_eQTL | ST22_GTEx_HighPPA_eQTL,
  ST31_STARNET_eQTL | ST28_STARNET_HighLD_eQTL | ST22_STARNET_HighPPA_eQTL
), na.rm = TRUE)]

gene_support[, AragamGeneSupportClass := fifelse(
  ST31_MostLikelyCausalGene,
  "Aragam prioritized causal gene",
  fifelse(
    ST31_PrimaryPoPS | ST31_SecondaryPoPS |
      ST31_GTEx_eQTL | ST31_STARNET_eQTL |
      ST27_GTEx_HighLD_eQTL | ST28_STARNET_HighLD_eQTL |
      ST22_GTEx_HighPPA_eQTL | ST22_STARNET_HighPPA_eQTL |
      ST31_ProteinAltering | ST31_DrugMREvidence,
    "Orthogonal CAD-specific gene support",
    fifelse(
      ST31_NearestGene,
      "Aragam nearest-gene support only",
      "No Aragam gene-table support"
    )
  )
)]

setorder(
  gene_support,
  -Any_AragamGeneSupport,
  -N_AragamGeneEvidenceAxes,
  Gene
)

# ============================================================
# 6. Gene-concordance at directly matched framework SNPs
# ============================================================

# This is the strongest V2G-oriented comparison in S6:
# for each Framework SNP + Framework gene combination, ask whether the same
# gene is supported by Aragam evidence attached to that SNP / CAD locus.

module_snp_gene <- rbindlist(
  lapply(seq_len(nrow(all_modules)), function(i) {
    snps <- parse_rsids(all_modules$SNPs[i])
    if (!length(snps)) return(NULL)
    data.table(
      FrameworkModuleID = all_modules$FrameworkModuleID[i],
      S5_Index = all_modules$S5_Index[i],
      Gene = all_modules$Gene[i],
      NatMed_DrugCategory = as.character(all_modules$NatMed_DrugCategory[i]),
      SNP = snps
    )
  }),
  fill = TRUE
)

# Build row-level maps required for exact gene concordance.
st27_gene_by_snp <- st27_long[, .(
  ST27_Genes_at_SNP = collapse_unique(ST27_Gene)
), by = SNP]

st28_gene_by_snp <- st28_long[, .(
  ST28_Genes_at_SNP = collapse_unique(ST28_Gene)
), by = SNP]

st22_gene_by_snp <- st22_long[, .(
  ST22_GTExGenes_at_SNP = collapse_unique(ST22_GTExGene),
  ST22_STARNETGenes_at_SNP = collapse_unique(ST22_STARNETGene)
), by = SNP]

st31_gene_by_snp <- st31_std[, .(
  ST31_MostLikelyGenes_at_SNP = collapse_unique(ST31_MostLikelyCausalGene),
  ST31_PrimaryPoPSGenes_at_SNP = collapse_unique(ST31_PrimaryPoPSGene),
  ST31_TopGTExGenes_at_SNP = collapse_unique(ST31_TopGTExGene),
  ST31_OtherGTExGenes_at_SNP = collapse_unique(ST31_OtherGTExGenes),
  ST31_TopSTARNETGenes_at_SNP = collapse_unique(ST31_TopSTARNETGene),
  ST31_OtherSTARNETGenes_at_SNP = collapse_unique(ST31_OtherSTARNETGenes),
  ST31_ProteinAlteringGenes_at_SNP = collapse_unique(ST31_ProteinAlteringGene),
  ST31_NPredictors_at_SNP = safe_max_int(ST31_NPredictors)
), by = SNP]

module_snp_gene <- merge(
  module_snp_gene,
  variant_support,
  by = "SNP",
  all.x = TRUE,
  sort = FALSE
)

for (lk in list(st27_gene_by_snp, st28_gene_by_snp, st22_gene_by_snp, st31_gene_by_snp)) {
  module_snp_gene <- merge(
    module_snp_gene,
    lk,
    by = "SNP",
    all.x = TRUE,
    sort = FALSE
  )
}

# Gene-concordance flags: exact gene-symbol token match only.
module_snp_gene[, GeneConcordant_ST27 := mapply(
  contains_gene, ST27_Genes_at_SNP, Gene
)]
module_snp_gene[, GeneConcordant_ST28 := mapply(
  contains_gene, ST28_Genes_at_SNP, Gene
)]
module_snp_gene[, GeneConcordant_ST22_GTEx := mapply(
  contains_gene, ST22_GTExGenes_at_SNP, Gene
)]
module_snp_gene[, GeneConcordant_ST22_STARNET := mapply(
  contains_gene, ST22_STARNETGenes_at_SNP, Gene
)]
module_snp_gene[, GeneConcordant_ST31_MostLikely := mapply(
  contains_gene, ST31_MostLikelyGenes_at_SNP, Gene
)]
module_snp_gene[, GeneConcordant_ST31_PoPS := mapply(
  contains_gene, ST31_PrimaryPoPSGenes_at_SNP, Gene
)]
module_snp_gene[, GeneConcordant_ST31_GTEx := mapply(
  function(a, b, g) contains_gene(c(a, b), g),
  ST31_TopGTExGenes_at_SNP,
  ST31_OtherGTExGenes_at_SNP,
  Gene
)]
module_snp_gene[, GeneConcordant_ST31_STARNET := mapply(
  function(a, b, g) contains_gene(c(a, b), g),
  ST31_TopSTARNETGenes_at_SNP,
  ST31_OtherSTARNETGenes_at_SNP,
  Gene
)]
module_snp_gene[, GeneConcordant_ST31_ProteinAltering := mapply(
  contains_gene, ST31_ProteinAlteringGenes_at_SNP, Gene
)]

module_snp_gene[, Any_GeneConcordant_AragamV2G :=
  GeneConcordant_ST27 |
    GeneConcordant_ST28 |
    GeneConcordant_ST22_GTEx |
    GeneConcordant_ST22_STARNET |
    GeneConcordant_ST31_MostLikely |
    GeneConcordant_ST31_PoPS |
    GeneConcordant_ST31_GTEx |
    GeneConcordant_ST31_STARNET |
    GeneConcordant_ST31_ProteinAltering]

setorder(
  module_snp_gene,
  -Any_GeneConcordant_AragamV2G,
  -Any_AragamVariantSupport,
  Gene,
  NatMed_DrugCategory,
  SNP
)

# ============================================================
# 7. Module-level post-prioritization summary
# ============================================================

module_agg <- module_snp_gene[, .(
  N_FrameworkSNPs = uniqueN(SNP),
  N_AragamMatchedSNPs = uniqueN(SNP[Any_AragamVariantSupport %in% TRUE]),
  N_FineMappedSNPs = uniqueN(SNP[Any_FineMappingSupport %in% TRUE]),
  N_ST21_CredibleSetSNPs = uniqueN(SNP[ST21_CredibleSet %in% TRUE]),
  N_ST22_HighPPATableSNPs = uniqueN(SNP[ST22_HighPPA_TableMatch %in% TRUE]),
  N_ST27_GTExLinkedSNPs = uniqueN(SNP[ST27_GTEx_HighLD_Support %in% TRUE]),
  N_ST28_STARNETLinkedSNPs = uniqueN(SNP[ST28_STARNET_HighLD_Support %in% TRUE]),
  N_ST31_CADSentinelSNPs = uniqueN(SNP[ST31_CADSentinelMatch %in% TRUE]),
  N_GeneConcordantAragamSNPs = uniqueN(SNP[Any_GeneConcordant_AragamV2G %in% TRUE]),
  Max_ST21_PPA = safe_max_num(ST21_MaxPPA),
  Max_ST22_PPA = safe_max_num(ST22_MaxPPA),
  Max_ST31_BestSNP_PPA = safe_max_num(ST31_MaxBestSNP_PPA),
  Max_ST31_NPredictors = safe_max_int(ST31_MaxNPredictors),
  AragamMatchedSNPs = collapse_unique(SNP[Any_AragamVariantSupport %in% TRUE]),
  FineMappedSNPs = collapse_unique(SNP[Any_FineMappingSupport %in% TRUE]),
  GeneConcordantAragamSNPs = collapse_unique(SNP[Any_GeneConcordant_AragamV2G %in% TRUE]),
  Aragam_ST31_MostLikelyGenes = collapse_unique(ST31_MostLikelyCausalGenes),
  Aragam_GTExGenes = collapse_unique(c(ST27_GTExGenes, ST22_GTExGenes, ST31_TopGTExGenes, ST31_OtherGTExGenes)),
  Aragam_STARNETGenes = collapse_unique(c(ST28_STARNETGenes, ST22_STARNETGenes, ST31_TopSTARNETGenes, ST31_OtherSTARNETGenes)),
  Aragam_CADSentinels = collapse_unique(c(ST27_CADSentinels, ST28_CADSentinels)),
  Any_DirectVariantSupport = any(Any_AragamVariantSupport %in% TRUE),
  Any_FineMappingSupport = any(Any_FineMappingSupport %in% TRUE),
  Any_HighLD_eQTLSupport = any(Any_HighLD_eQTLSupport %in% TRUE),
  Any_GeneConcordant_AragamV2G = any(Any_GeneConcordant_AragamV2G %in% TRUE)
), by = .(
  FrameworkModuleID,
  S5_Index,
  Gene,
  NatMed_DrugCategory
)]

module_summary <- merge(
  all_modules,
  module_agg,
  by = c("FrameworkModuleID", "S5_Index", "Gene", "NatMed_DrugCategory"),
  all.x = TRUE,
  sort = FALSE
)

module_summary <- merge(
  module_summary,
  gene_support,
  by = "Gene",
  all.x = TRUE,
  sort = FALSE
)

# Modules without parsed rsIDs (unexpected, but handled safely).
for (nm in c(
  "Any_DirectVariantSupport",
  "Any_FineMappingSupport",
  "Any_HighLD_eQTLSupport",
  "Any_GeneConcordant_AragamV2G",
  "Any_AragamGeneSupport"
)) {
  if (nm %in% names(module_summary)) {
    module_summary[is.na(get(nm)), (nm) := FALSE]
  }
}

module_summary[, AragamCADSupportClass := fifelse(
  Any_GeneConcordant_AragamV2G,
  "Locus-reinforced: direct variant/locus + gene-concordant CAD support",
  fifelse(
    Any_DirectVariantSupport & Any_AragamGeneSupport,
    "Direct variant/locus support + separate gene-level CAD support",
    fifelse(
      Any_DirectVariantSupport,
      "Direct variant/locus CAD support; assigned gene not concordant in matched Aragam evidence",
      fifelse(
        Any_AragamGeneSupport,
        "Gene-level CAD support only; no direct framework-SNP match in Aragam tables",
        "No Aragam support identified in queried supplementary resources"
      )
    )
  )
)]

# Descriptive-only evidence-axis count; not used for ranking or tier assignment.
module_summary[, N_AragamModuleSupportAxes := rowSums(cbind(
  Any_FineMappingSupport,
  N_ST27_GTExLinkedSNPs > 0,
  N_ST28_STARNETLinkedSNPs > 0,
  N_ST31_CADSentinelSNPs > 0,
  Any_GeneConcordant_AragamV2G,
  Any_AragamGeneSupport
), na.rm = TRUE)]

# Keep primary S5 fields first, then Aragam support fields.
front_cols <- c(
  "S5_Index",
  "FrameworkModuleID",
  "Gene",
  "NatMed_DrugCategory",
  "NatMed_CategoryCode",
  "CandidateClass",
  "CandidateTierA",
  "CandidateTierB",
  "CandidateTierC",
  "UniqueSNPs",
  "Tier1SNPs",
  "UniqueMetabolites",
  "Tier1Metabolites",
  "SNPs",
  "Metabolites",
  "N_FrameworkSNPs",
  "N_AragamMatchedSNPs",
  "N_FineMappedSNPs",
  "N_ST21_CredibleSetSNPs",
  "N_ST22_HighPPATableSNPs",
  "N_ST27_GTExLinkedSNPs",
  "N_ST28_STARNETLinkedSNPs",
  "N_ST31_CADSentinelSNPs",
  "N_GeneConcordantAragamSNPs",
  "AragamMatchedSNPs",
  "FineMappedSNPs",
  "GeneConcordantAragamSNPs",
  "Max_ST21_PPA",
  "Max_ST22_PPA",
  "Max_ST31_BestSNP_PPA",
  "Max_ST31_NPredictors",
  "Aragam_ST31_MostLikelyGenes",
  "Aragam_GTExGenes",
  "Aragam_STARNETGenes",
  "Aragam_CADSentinels",
  "Any_DirectVariantSupport",
  "Any_FineMappingSupport",
  "Any_HighLD_eQTLSupport",
  "Any_GeneConcordant_AragamV2G",
  "Any_AragamGeneSupport",
  "AragamGeneSupportClass",
  "N_AragamModuleSupportAxes",
  "AragamCADSupportClass"
)
front_cols <- intersect(front_cols, names(module_summary))
rest_cols <- setdiff(names(module_summary), front_cols)
setcolorder(module_summary, c(front_cols, rest_cols))

# ============================================================
# 8. Tier-specific annotated tables
# ============================================================

# Match by framework Gene + drug category. These are the canonical module keys.
annotate_tier <- function(tier_dt) {
  z <- copy(tier_dt)
  z[, Gene := norm_gene(Gene)]

  keep_ann <- module_summary[, .(
    Gene,
    NatMed_DrugCategory,
    S5_Index,
    FrameworkModuleID,
    N_FrameworkSNPs,
    N_AragamMatchedSNPs,
    N_FineMappedSNPs,
    N_ST21_CredibleSetSNPs,
    N_ST22_HighPPATableSNPs,
    N_ST27_GTExLinkedSNPs,
    N_ST28_STARNETLinkedSNPs,
    N_ST31_CADSentinelSNPs,
    N_GeneConcordantAragamSNPs,
    AragamMatchedSNPs,
    FineMappedSNPs,
    GeneConcordantAragamSNPs,
    Max_ST21_PPA,
    Max_ST22_PPA,
    Max_ST31_BestSNP_PPA,
    Max_ST31_NPredictors,
    Aragam_ST31_MostLikelyGenes,
    Aragam_GTExGenes,
    Aragam_STARNETGenes,
    Aragam_CADSentinels,
    Any_DirectVariantSupport,
    Any_FineMappingSupport,
    Any_HighLD_eQTLSupport,
    Any_GeneConcordant_AragamV2G,
    Any_AragamGeneSupport,
    AragamGeneSupportClass,
    N_AragamModuleSupportAxes,
    AragamCADSupportClass
  )]

  merge(
    z,
    keep_ann,
    by = c("Gene", "NatMed_DrugCategory"),
    all.x = TRUE,
    sort = FALSE
  )
}

TierA <- annotate_tier(tierA_s5)
TierB <- annotate_tier(tierB_s5)
TierC <- annotate_tier(tierC_s5)

# ============================================================
# 9. QC summary
# ============================================================

qc <- data.table(
  Metric = c(
    "S5 all modules",
    "S5 unique genes",
    "S5 unique framework SNPs",
    "Framework SNPs matching any queried Aragam variant table",
    "Framework SNPs in ST21 95% credible sets",
    "Framework SNPs matching ST22 high-PPA table",
    "Framework SNPs matching ST27 CAD-sentinel/high-LD GTEx table",
    "Framework SNPs matching ST28 CAD-sentinel/high-LD STARNET table",
    "Framework SNPs matching ST31 CAD sentinel",
    "Framework genes with any Aragam gene-level support",
    "Modules with direct Aragam variant/locus support",
    "Modules with fine-mapping support",
    "Modules with gene-concordant Aragam V2G support at matched SNP/locus",
    "Modules with Aragam gene-level support",
    "Tier A modules",
    "Tier A with direct Aragam variant/locus support",
    "Tier A with gene-concordant Aragam V2G support",
    "Tier A with any Aragam gene-level support"
  ),
  Value = c(
    nrow(all_modules),
    length(framework_genes),
    length(all_framework_snps),
    variant_support[Any_AragamVariantSupport == TRUE, uniqueN(SNP)],
    variant_support[ST21_CredibleSet == TRUE, uniqueN(SNP)],
    variant_support[ST22_HighPPA_TableMatch == TRUE, uniqueN(SNP)],
    variant_support[ST27_GTEx_HighLD_Support == TRUE, uniqueN(SNP)],
    variant_support[ST28_STARNET_HighLD_Support == TRUE, uniqueN(SNP)],
    variant_support[ST31_CADSentinelMatch == TRUE, uniqueN(SNP)],
    gene_support[Any_AragamGeneSupport == TRUE, uniqueN(Gene)],
    module_summary[Any_DirectVariantSupport == TRUE, .N],
    module_summary[Any_FineMappingSupport == TRUE, .N],
    module_summary[Any_GeneConcordant_AragamV2G == TRUE, .N],
    module_summary[Any_AragamGeneSupport == TRUE, .N],
    nrow(TierA),
    TierA[Any_DirectVariantSupport == TRUE, .N],
    TierA[Any_GeneConcordant_AragamV2G == TRUE, .N],
    TierA[Any_AragamGeneSupport == TRUE, .N]
  )
)

fwrite(qc, QC_TSV, sep = "\t", quote = FALSE, na = "NA")

# ============================================================
# 10. Dictionary / interpretation guardrails
# ============================================================

dictionary <- data.table(
  Sheet = c(
    "ModuleSummary",
    "TierA / TierB / TierC",
    "VariantSupport",
    "ModuleVariantV2G",
    "GeneSupport",
    "QC"
  ),
  Unit = c(
    "Framework Gene + NatMed drug-category module",
    "Framework candidate module",
    "Unique Framework mGWAS instrument rsID",
    "Framework module + instrument rsID",
    "Unique Framework-assigned gene",
    "Summary metric"
  ),
  Description = c(
    paste0(
      "All S5 modules annotated post-prioritization with Aragam CAD-specific ",
      "fine-mapping, high-LD eQTL and causal-gene evidence. Original candidate tiers are unchanged."
    ),
    paste0(
      "Original S5 Tier A/B/C rows plus the same Aragam annotations used in ModuleSummary. ",
      "Aragam support is not used to create or re-rank tiers."
    ),
    paste0(
      "One row per unique Framework SNP. Direct rsID overlap is tracked against ST21, ST22, ",
      "ST27, ST28 and ST31. ST27/ST28 indicate high-LD eQTL support, not formal eQTL colocalization."
    ),
    paste0(
      "One row per Framework module-SNP pair. Compares the VEP-assigned Framework gene with ",
      "Aragam genes supported at the directly matched SNP/CAD locus."
    ),
    paste0(
      "Gene-level CAD support irrespective of direct Framework-SNP overlap. This is kept separate ",
      "from locus-specific support to avoid conflating network-level propagation with CAD-locus anchoring."
    ),
    "Counts used to audit the integration."
  )
)

field_dictionary <- data.table(
  Field = c(
    "Any_DirectVariantSupport",
    "Any_FineMappingSupport",
    "Any_HighLD_eQTLSupport",
    "Any_GeneConcordant_AragamV2G",
    "Any_AragamGeneSupport",
    "N_AragamSupportAxes",
    "N_AragamModuleSupportAxes",
    "AragamCADSupportClass",
    "ST21_CredibleSet",
    "ST22_HighPPA_TableMatch",
    "ST27_GTEx_HighLD_Support",
    "ST28_STARNET_HighLD_Support",
    "ST31_CADSentinelMatch"
  ),
  Meaning = c(
    "Framework SNP appears directly in at least one queried Aragam variant/locus table.",
    "Framework SNP appears in ST21 95% credible-set variants and/or ST22 high-PPA table.",
    "Framework SNP appears as a CAD sentinel or high-LD eQTL proxy/index SNP in ST27 or ST28.",
    "At a directly matched Framework SNP/CAD locus, Aragam evidence supports the same gene assigned by Framework VEP.",
    "Framework-assigned gene appears anywhere in queried Aragam CAD gene-prioritization/eQTL resources, regardless of direct SNP match.",
    "Descriptive count of Aragam variant-level support axes. NOT a prioritization score.",
    "Descriptive count of module-level Aragam evidence axes. NOT a prioritization score.",
    "Textual evidence class distinguishing locus-reinforced, variant-supported and gene-only support.",
    "Direct rsID match to an ST21 95% credible-set variant.",
    "Direct rsID match to an ST22 high-PPA/fine-mapping-eQTL table variant field.",
    "Direct rsID match to ST27 CAD sentinel or GTEx high-LD eQTL proxy.",
    "Direct rsID match to ST28 CAD sentinel or STARNET high-LD index eSNP.",
    "Direct rsID match to the CAD association sentinel in ST31."
  ),
  InterpretationGuardrail = c(
    "Direct table overlap does not by itself establish causal drug action.",
    "Fine-mapping support strengthens CAD-locus anchoring but is not required by the core framework.",
    "Do not describe ST27/ST28 as formal eQTL colocalization.",
    "Strongest V2G-oriented concordance in S6, but still post-prioritization support.",
    "Gene-level support is weaker than same-locus gene-concordant evidence.",
    "Do not use as a biological ranking score.",
    "Do not use as a biological ranking score.",
    "Evidence class, not a causal or clinical tier.",
    "Absence from ST21 does not mean the variant is not functional or disease-relevant.",
    "ST22 is a selected high-PPA subset and should not be treated as genome-wide coverage.",
    "High LD with a CAD sentinel is supportive but not equivalent to fine-mapping or colocalization.",
    "High LD with a CAD sentinel is supportive but not equivalent to fine-mapping or colocalization.",
    "Sentinel identity is direct locus overlap, not proof that the assigned gene is causal."
  )
)

# ============================================================
# 11. Write Supplementary_Data_S6.xlsx
# ============================================================

wb <- openxlsx::createWorkbook(creator = "Framework_v1")

header_style <- openxlsx::createStyle(
  fontColour = "#FFFFFF",
  fgFill = "#2F5597",
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  wrapText = TRUE,
  border = "Bottom",
  borderColour = "#D9E2F3"
)

subheader_style <- openxlsx::createStyle(
  fontColour = "#FFFFFF",
  fgFill = "#5B9BD5",
  textDecoration = "bold",
  halign = "left",
  valign = "center",
  wrapText = TRUE
)

wrap_style <- openxlsx::createStyle(
  valign = "top",
  wrapText = TRUE
)

bool_true_style <- openxlsx::createStyle(
  fgFill = "#E2F0D9",
  fontColour = "#375623"
)

fine_map_style <- openxlsx::createStyle(
  fgFill = "#FFF2CC",
  fontColour = "#7F6000"
)

locus_reinforced_style <- openxlsx::createStyle(
  fgFill = "#DDEBF7",
  fontColour = "#1F4E78"
)

add_table_sheet <- function(wb, name, dt, freeze_first_row = TRUE) {
  openxlsx::addWorksheet(wb, name)
  openxlsx::writeData(
    wb,
    sheet = name,
    x = as.data.frame(dt),
    headerStyle = header_style,
    withFilter = TRUE,
    keepNA = TRUE,
    na.string = "NA"
  )

  if (freeze_first_row) openxlsx::freezePane(wb, name, firstRow = TRUE)

  openxlsx::addStyle(
    wb,
    sheet = name,
    style = wrap_style,
    rows = seq_len(nrow(dt) + 1L),
    cols = seq_len(ncol(dt)),
    gridExpand = TRUE,
    stack = TRUE
  )

  # Conservative widths: avoid enormous columns from long SNP/gene lists.
  if (ncol(dt) > 0L) {
    widths <- rep(14, ncol(dt))
    nm <- names(dt)
    widths[grepl("Gene$|^Gene$|SNP$|^SNP$|Index|Count|^N_", nm)] <- 13
    widths[grepl("PPA|r2|R2|FDR|P$|Predictor", nm)] <- 14
    widths[grepl("Category|Class|Support|Tissue|Loci|Genes|SNPs|Metabolites|Role|Reason", nm)] <- 26
    widths[grepl("Description|Meaning|Guardrail", nm)] <- 42
    widths <- pmin(widths, 42)
    openxlsx::setColWidths(wb, name, cols = seq_along(widths), widths = widths)
  }

  # Conditional formatting where the fields exist.
  true_cols <- which(names(dt) %in% c(
    "Any_DirectVariantSupport",
    "Any_FineMappingSupport",
    "Any_HighLD_eQTLSupport",
    "Any_GeneConcordant_AragamV2G",
    "Any_AragamGeneSupport",
    "ST21_CredibleSet",
    "ST22_HighPPA_TableMatch",
    "ST27_GTEx_HighLD_Support",
    "ST28_STARNET_HighLD_Support",
    "ST31_CADSentinelMatch"
  ))

  for (cc in true_cols) {
    if (nrow(dt) > 0L) {
      openxlsx::conditionalFormatting(
        wb,
        name,
        cols = cc,
        rows = 2:(nrow(dt) + 1L),
        rule = paste0(openxlsx::int2col(cc), "2=TRUE"),
        style = bool_true_style,
        type = "expression"
      )
    }
  }

  support_class_col <- which(names(dt) == "AragamCADSupportClass")
  if (length(support_class_col) && nrow(dt) > 0L) {
    cc <- support_class_col[1]
    col_letter <- openxlsx::int2col(cc)
    openxlsx::conditionalFormatting(
      wb, name, cols = cc, rows = 2:(nrow(dt) + 1L),
      rule = paste0('ISNUMBER(SEARCH("Locus-reinforced",', col_letter, '2))'),
      style = locus_reinforced_style,
      type = "expression"
    )
  }

  fine_col <- which(names(dt) == "Any_FineMappingSupport")
  if (length(fine_col) && nrow(dt) > 0L) {
    cc <- fine_col[1]
    openxlsx::conditionalFormatting(
      wb, name, cols = cc, rows = 2:(nrow(dt) + 1L),
      rule = paste0(openxlsx::int2col(cc), "2=TRUE"),
      style = fine_map_style,
      type = "expression"
    )
  }
}

# Dictionary sheet with two blocks.
openxlsx::addWorksheet(wb, "Dictionary")
openxlsx::writeData(wb, "Dictionary", "Supplementary Data S6 — Aragam CAD-specific post-prioritization support", startRow = 1, startCol = 1)
openxlsx::mergeCells(wb, "Dictionary", cols = 1:3, rows = 1)
openxlsx::addStyle(wb, "Dictionary", subheader_style, rows = 1, cols = 1:3, gridExpand = TRUE)

openxlsx::writeData(wb, "Dictionary", dictionary, startRow = 3, startCol = 1, headerStyle = header_style)
field_start <- 3 + nrow(dictionary) + 3
openxlsx::writeData(wb, "Dictionary", "Key derived fields", startRow = field_start - 1, startCol = 1)
openxlsx::mergeCells(wb, "Dictionary", cols = 1:3, rows = field_start - 1)
openxlsx::addStyle(wb, "Dictionary", subheader_style, rows = field_start - 1, cols = 1:3, gridExpand = TRUE)
openxlsx::writeData(wb, "Dictionary", field_dictionary, startRow = field_start, startCol = 1, headerStyle = header_style)
openxlsx::setColWidths(wb, "Dictionary", cols = 1, widths = 30)
openxlsx::setColWidths(wb, "Dictionary", cols = 2:3, widths = 55)
openxlsx::addStyle(
  wb, "Dictionary", wrap_style,
  rows = 1:(field_start + nrow(field_dictionary) + 2), cols = 1:3,
  gridExpand = TRUE, stack = TRUE
)
openxlsx::freezePane(wb, "Dictionary", firstRow = TRUE)

add_table_sheet(wb, "ModuleSummary", module_summary)
add_table_sheet(wb, "TierA", TierA)
add_table_sheet(wb, "TierB", TierB)
add_table_sheet(wb, "TierC", TierC)
add_table_sheet(wb, "VariantSupport", variant_support)
add_table_sheet(wb, "ModuleVariantV2G", module_snp_gene)
add_table_sheet(wb, "GeneSupport", gene_support)
add_table_sheet(wb, "QC", qc)

# Number formats for PPA/r2-like columns.
for (sheet_name in c("ModuleSummary", "TierA", "TierB", "TierC", "VariantSupport", "ModuleVariantV2G")) {
  dt <- switch(
    sheet_name,
    ModuleSummary = module_summary,
    TierA = TierA,
    TierB = TierB,
    TierC = TierC,
    VariantSupport = variant_support,
    ModuleVariantV2G = module_snp_gene
  )

  num_cols <- which(grepl("PPA|R2|r2", names(dt)))
  if (length(num_cols) && nrow(dt) > 0L) {
    openxlsx::addStyle(
      wb,
      sheet = sheet_name,
      style = openxlsx::createStyle(numFmt = "0.000"),
      rows = 2:(nrow(dt) + 1L),
      cols = num_cols,
      gridExpand = TRUE,
      stack = TRUE
    )
  }
}

# Metadata sheet keeps provenance explicit.
metadata <- data.table(
  Field = c(
    "Framework input",
    "Aragam input",
    "Aragam tables used",
    "Output purpose",
    "Candidate-tier effect",
    "Important terminology",
    "Generated on"
  ),
  Value = c(
    S5_FILE,
    ARAGAM_FILE,
    "ST20; ST21; ST22; ST27; ST28; ST31",
    "Post-prioritization CAD-specific fine-mapping / V2G / high-LD eQTL annotation",
    "None — Tier A/B/C are inherited from Supplementary_Data_S5.xlsx and are not re-ranked",
    "ST27/ST28 = high-LD eQTL support, not formal eQTL colocalization",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )
)
add_table_sheet(wb, "Metadata", metadata)

openxlsx::saveWorkbook(
  wb,
  OUTPUT_FILE,
  overwrite = TRUE
)

# ============================================================
# 12. Final console report
# ============================================================

cat("\n============================================================\n")
cat("Supplementary_Data_S6.xlsx created successfully\n")
cat("============================================================\n\n")
print(qc)

cat("\nWorkbook sheets:\n")
cat("  Dictionary\n")
cat("  ModuleSummary\n")
cat("  TierA\n")
cat("  TierB\n")
cat("  TierC\n")
cat("  VariantSupport\n")
cat("  ModuleVariantV2G\n")
cat("  GeneSupport\n")
cat("  QC\n")
cat("  Metadata\n")

cat("\nOutput workbook:\n  ", OUTPUT_FILE, "\n", sep = "")
cat("QC TSV:\n  ", QC_TSV, "\n", sep = "")

cat("\nInterpretation:\n")
cat("  Aragam evidence is an orthogonal CAD-specific post-prioritization\n")
cat("  annotation layer. It does not alter Framework_v1 candidate tiers.\n")
cat("  Direct SNP/locus support and gene-level support are deliberately\n")
cat("  reported separately.\n")
