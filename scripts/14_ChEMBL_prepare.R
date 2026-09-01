#!/usr/bin/env Rscript

# ============================================================
# Framework_v1
# 14_ChEMBL_prepare.R
#
# Download and freeze a ChEMBL mechanism-of-action snapshot for
# the already accepted DGIdb -> NatMed drug mappings.
#
# Inputs:
#   data/processed/DGIdb_to_NatMed_DrugCategory_Mapping.tsv
#
# Outputs:
#   data/raw/ChEMBL/ChEMBL_drug_mapping.tsv
#   data/raw/ChEMBL/ChEMBL_drug_mechanisms.tsv
#   data/raw/ChEMBL/ChEMBL_targets.tsv
#   data/raw/ChEMBL/ChEMBL_target_components.tsv
#   data/raw/ChEMBL/ChEMBL_download_metadata.tsv
#   data/processed/ChEMBL_MoA_GeneDrugEvidence.tsv
#   data/processed/ChEMBL_Preparation_QC.tsv
#
# This script prepares data only. It does NOT benchmark our
# gene-drug modules yet.
# ============================================================

options(stringsAsFactors = FALSE, warn = 1, scipen = 999)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
source(file.path(dirname(script_path), "..", "R", "bootstrap.R"))

for (pkg in c("data.table", "jsonlite", "httr")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required R package missing: ", pkg, call. = FALSE)
  }
}
library(data.table)

BASE <- "https://www.ebi.ac.uk/chembl/api/data"
raw_dir <- file.path(DIR_RAW, "ChEMBL")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

mapping_file <- file.path(
  DIR_PROCESSED, "DGIdb_to_NatMed_DrugCategory_Mapping.tsv"
)

if (!file.exists(mapping_file)) {
  stop("Missing mapping file: ", mapping_file, call. = FALSE)
}

norm_drug <- function(x) {
  x <- tolower(trimws(as.character(x)))
  gsub("[[:space:]]+", " ", x)
}

api_get <- function(url, query = list(), attempts = 4L) {
  for (i in seq_len(attempts)) {
    r <- try(
      httr::GET(
        url,
        query = query,
        httr::timeout(60),
        httr::user_agent("Framework_v1 ChEMBL validation")
      ),
      silent = TRUE
    )
    if (!inherits(r, "try-error") && httr::status_code(r) == 200L) {
      txt <- httr::content(r, as = "text", encoding = "UTF-8")
      return(jsonlite::fromJSON(txt, simplifyVector = FALSE))
    }
    Sys.sleep(2 * i)
  }
  stop("ChEMBL API request failed: ", url, call. = FALSE)
}

get_all <- function(endpoint, query = list(), key) {
  query$format <- "json"
  query$limit <- 1000
  query$offset <- 0
  out <- list()

  repeat {
    obj <- api_get(paste0(BASE, "/", endpoint, ".json"), query)
    rows <- obj[[key]]
    if (length(rows)) out <- c(out, rows)

    meta <- obj$page_meta
    next_page <- meta[["next"]]
    if (is.null(next_page) || length(next_page) == 0L ||
        is.na(next_page) || identical(next_page, "")) break
    query$offset <- query$offset + query$limit
  }
  out
}

scalar <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA_character_)
  as.character(x[[1]])
}

collapse_field <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA_character_)
  vals <- unique(unlist(x, recursive = TRUE, use.names = FALSE))
  vals <- vals[!is.na(vals) & vals != ""]
  if (!length(vals)) NA_character_ else paste(vals, collapse = ";")
}

cat("\n=== ChEMBL mechanism-of-action preparation ===\n")

mp <- fread(mapping_file, sep = "\t", na.strings = c("", "NA"))
need <- c(
  "DGIdb_Drug", "DGIdb_Drug_Normalized", "MappingStatus",
  "NatMed_DrugCategory", "NatMed_CategoryCode"
)
miss <- setdiff(need, names(mp))
if (length(miss)) stop("Mapping missing columns: ", paste(miss, collapse=", "))

accepted <- unique(mp[
  MappingStatus %in% c("REUSED_CPIC", "MANUAL_ACCEPTED") &
    !is.na(NatMed_DrugCategory),
  .(
    InputDrug = DGIdb_Drug,
    InputDrugNormalized = norm_drug(DGIdb_Drug_Normalized),
    NatMed_DrugCategory,
    NatMed_CategoryCode,
    MappingStatus
  )
])

cat("Accepted mapped drugs: ", nrow(accepted), "\n", sep="")

# ------------------------------------------------------------
# Detect ChEMBL release from API
# ------------------------------------------------------------
release_obj <- api_get(paste0(BASE, "/chembl_release.json"),
                       list(format="json", limit=1000))
release_rows <- release_obj$chembl_releases
release_dt <- rbindlist(lapply(release_rows, function(z) {
  data.table(
    chembl_release = scalar(z$chembl_release),
    creation_date = scalar(z$creation_date)
  )
}), fill=TRUE)

release_dt[, relnum := suppressWarnings(
  as.integer(gsub("[^0-9]", "", chembl_release))
)]
setorder(release_dt, -relnum)
current_release <- release_dt[1]

cat("ChEMBL API release: ", current_release$chembl_release, "\n", sep="")

# ------------------------------------------------------------
# Drug -> ChEMBL molecule mapping
# Search by name, then require exact normalized pref_name or synonym.
# ------------------------------------------------------------
drug_map_rows <- list()

for (i in seq_len(nrow(accepted))) {
  d <- accepted[i]
  cat(sprintf("[%d/%d] molecule search: %s\n",
              i, nrow(accepted), d$InputDrug))

  res <- get_all(
    "molecule/search",
    query = list(q = d$InputDrug),
    key = "molecules"
  )

  candidates <- list()

  if (length(res)) {
    for (m in res) {
      chembl_id <- scalar(m$molecule_chembl_id)
      pref <- scalar(m$pref_name)

      syns <- character()
      if (!is.null(m$molecule_synonyms) && length(m$molecule_synonyms)) {
        syns <- unique(unlist(lapply(m$molecule_synonyms, function(s) {
          c(s$molecule_synonym, s$synonyms)
        }), use.names = FALSE))
      }

      names_norm <- unique(norm_drug(c(pref, syns)))
      exact <- d$InputDrugNormalized %in% names_norm

      candidates[[length(candidates)+1L]] <- data.table(
        InputDrug = d$InputDrug,
        InputDrugNormalized = d$InputDrugNormalized,
        NatMed_DrugCategory = d$NatMed_DrugCategory,
        NatMed_CategoryCode = d$NatMed_CategoryCode,
        MappingStatus = d$MappingStatus,
        MoleculeChEMBLID = chembl_id,
        ChEMBLPrefName = pref,
        ExactNameOrSynonym = exact,
        MoleculeType = scalar(m$molecule_type),
        MaxPhase = scalar(m$max_phase)
      )
    }
  }

  cand <- if (length(candidates)) rbindlist(candidates, fill=TRUE) else data.table()

  if (nrow(cand) && any(cand$ExactNameOrSynonym %in% TRUE)) {
    cand <- cand[ExactNameOrSynonym == TRUE]
    # Prefer higher max phase if multiple exact salt/parent candidates.
    cand[, phase_num := suppressWarnings(as.numeric(MaxPhase))]
    setorder(cand, -phase_num, MoleculeChEMBLID)
    cand <- cand[1]
    cand[, phase_num := NULL]
    cand[, ChEMBLMappingStatus := "EXACT_NAME_OR_SYNONYM"]
  } else {
    cand <- data.table(
      InputDrug = d$InputDrug,
      InputDrugNormalized = d$InputDrugNormalized,
      NatMed_DrugCategory = d$NatMed_DrugCategory,
      NatMed_CategoryCode = d$NatMed_CategoryCode,
      MappingStatus = d$MappingStatus,
      MoleculeChEMBLID = NA_character_,
      ChEMBLPrefName = NA_character_,
      ExactNameOrSynonym = FALSE,
      MoleculeType = NA_character_,
      MaxPhase = NA_character_,
      ChEMBLMappingStatus = "NOT_EXACTLY_MAPPED"
    )
  }

  drug_map_rows[[i]] <- cand
  Sys.sleep(0.05)
}

drug_map <- rbindlist(drug_map_rows, fill=TRUE)
fwrite(drug_map,
       file.path(raw_dir, "ChEMBL_drug_mapping.tsv"),
       sep="\t", quote=FALSE, na="NA")

mapped_ids <- unique(na.omit(drug_map$MoleculeChEMBLID))
cat("Exactly mapped ChEMBL molecules: ", length(mapped_ids), "\n", sep="")

# ------------------------------------------------------------
# Molecule -> mechanism
# ------------------------------------------------------------
mech_rows <- list()

for (i in seq_along(mapped_ids)) {
  id <- mapped_ids[i]
  cat(sprintf("[%d/%d] mechanism: %s\n", i, length(mapped_ids), id))

  res <- get_all(
    "mechanism",
    query = list(molecule_chembl_id = id),
    key = "mechanisms"
  )

  if (length(res)) {
    mech_rows[[length(mech_rows)+1L]] <- rbindlist(
      lapply(res, function(z) {
        data.table(
          MoleculeChEMBLID = scalar(z$molecule_chembl_id),
          TargetChEMBLID = scalar(z$target_chembl_id),
          MechanismOfAction = scalar(z$mechanism_of_action),
          ActionType = scalar(z$action_type),
          DirectInteraction = scalar(z$direct_interaction),
          MolecularMechanism = scalar(z$molecular_mechanism),
          DiseaseEfficacy = scalar(z$disease_efficacy),
          MechanismComment = scalar(z$mechanism_comment),
          MechanismRefs = collapse_field(z$mechanism_refs)
        )
      }),
      fill=TRUE
    )
  }
  Sys.sleep(0.05)
}

mechanisms <- if (length(mech_rows)) rbindlist(mech_rows, fill=TRUE) else data.table()
fwrite(mechanisms,
       file.path(raw_dir, "ChEMBL_drug_mechanisms.tsv"),
       sep="\t", quote=FALSE, na="NA")

# ------------------------------------------------------------
# Targets and components
# ChEMBL target objects contain target_components; preserve them.
# ------------------------------------------------------------
target_ids <- if (nrow(mechanisms)) unique(na.omit(mechanisms$TargetChEMBLID)) else character()
target_rows <- list()
component_rows <- list()

for (i in seq_along(target_ids)) {
  tid <- target_ids[i]
  cat(sprintf("[%d/%d] target: %s\n", i, length(target_ids), tid))

  z <- api_get(paste0(BASE, "/target/", tid, ".json"),
               list(format="json"))

  target_rows[[i]] <- data.table(
    TargetChEMBLID = scalar(z$target_chembl_id),
    TargetPrefName = scalar(z$pref_name),
    TargetType = scalar(z$target_type),
    Organism = scalar(z$organism),
    TaxID = scalar(z$tax_id)
  )

  comps <- z$target_components
  if (!is.null(comps) && length(comps)) {
    for (c in comps) {
      acc <- scalar(c$accession)
      cid <- scalar(c$component_id)

      syns <- c$target_component_synonyms
      gene_symbols <- character()
      all_synonyms <- character()

      if (!is.null(syns) && length(syns)) {
        for (s in syns) {
          sval <- scalar(s$component_synonym)
          stype <- toupper(scalar(s$syn_type))
          if (!is.na(sval)) all_synonyms <- c(all_synonyms, sval)
          if (!is.na(stype) &&
              stype %in% c("GENE_SYMBOL", "HGNC_SYMBOL")) {
            gene_symbols <- c(gene_symbols, sval)
          }
        }
      }

      component_rows[[length(component_rows)+1L]] <- data.table(
        TargetChEMBLID = tid,
        ComponentID = cid,
        Accession = acc,
        GeneSymbols = if (length(gene_symbols))
          paste(sort(unique(gene_symbols)), collapse=";") else NA_character_,
        ComponentSynonyms = if (length(all_synonyms))
          paste(sort(unique(all_synonyms)), collapse=";") else NA_character_
      )
    }
  }
  Sys.sleep(0.05)
}

targets <- if (length(target_rows)) rbindlist(target_rows, fill=TRUE) else data.table()
components <- if (length(component_rows)) rbindlist(component_rows, fill=TRUE) else data.table()

fwrite(targets,
       file.path(raw_dir, "ChEMBL_targets.tsv"),
       sep="\t", quote=FALSE, na="NA")
fwrite(components,
       file.path(raw_dir, "ChEMBL_target_components.tsv"),
       sep="\t", quote=FALSE, na="NA")

# ------------------------------------------------------------
# Process human target-component genes
# ------------------------------------------------------------
evidence <- data.table()

if (nrow(mechanisms) && nrow(targets) && nrow(components)) {
  x <- merge(mechanisms, targets, by="TargetChEMBLID", all.x=TRUE)
  x <- x[Organism == "Homo sapiens"]

  x <- merge(x, components, by="TargetChEMBLID", all.x=TRUE)
  x <- merge(x, drug_map[
    !is.na(MoleculeChEMBLID),
    .(
      MoleculeChEMBLID,
      InputDrug,
      NatMed_DrugCategory,
      NatMed_CategoryCode
    )
  ], by="MoleculeChEMBLID", all.x=TRUE)

  # Expand semicolon-delimited HGNC/gene symbols.
  ev <- x[!is.na(GeneSymbols) & GeneSymbols != ""]
  if (nrow(ev)) {
    ev <- ev[, {
      gs <- trimws(unlist(strsplit(GeneSymbols, ";", fixed=TRUE)))
      .(Gene = toupper(gs))
    }, by = setdiff(names(ev), "GeneSymbols")]

    evidence <- unique(ev[, .(
      Gene,
      InputDrug,
      NatMed_DrugCategory,
      NatMed_CategoryCode,
      MoleculeChEMBLID,
      ChEMBLPrefName = NA_character_,
      TargetChEMBLID,
      TargetPrefName,
      TargetType,
      ActionType,
      MechanismOfAction,
      DirectInteraction,
      MolecularMechanism,
      DiseaseEfficacy,
      Accession,
      MechanismRefs
    )])
  }
}

fwrite(
  evidence,
  file.path(DIR_PROCESSED, "ChEMBL_MoA_GeneDrugEvidence.tsv"),
  sep="\t", quote=FALSE, na="NA"
)

# ------------------------------------------------------------
# Metadata and QC
# ------------------------------------------------------------
metadata <- data.table(
  Field = c(
    "API_base",
    "ChEMBL_release",
    "ChEMBL_release_creation_date",
    "download_date",
    "mapping_input"
  ),
  Value = c(
    BASE,
    current_release$chembl_release,
    current_release$creation_date,
    format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    mapping_file
  )
)

fwrite(metadata,
       file.path(raw_dir, "ChEMBL_download_metadata.tsv"),
       sep="\t", quote=FALSE)

n_mapped_drugs <- uniqueN(drug_map[!is.na(MoleculeChEMBLID), InputDrug])
drugs_with_moa <- if (nrow(mechanisms)) uniqueN(
  mechanisms$MoleculeChEMBLID
) else 0L
human_targets <- if (nrow(targets)) uniqueN(
  targets[Organism == "Homo sapiens", TargetChEMBLID]
) else 0L
human_genes <- if (nrow(evidence)) uniqueN(evidence$Gene) else 0L
gene_drug_cat <- if (nrow(evidence)) uniqueN(
  paste(evidence$Gene, evidence$NatMed_DrugCategory, sep="||")
) else 0L

qc <- data.table(
  Metric = c(
    "Accepted input drugs",
    "Input drugs exactly mapped to ChEMBL",
    "Input drugs not exactly mapped",
    "Unique ChEMBL molecules mapped",
    "Mapped ChEMBL molecules with >=1 mechanism",
    "Mechanism rows downloaded",
    "Unique mechanism targets",
    "Human mechanism targets",
    "Human target-component genes",
    "Human Gene-NatMed category MoA units"
  ),
  Value = c(
    nrow(accepted),
    n_mapped_drugs,
    nrow(accepted) - n_mapped_drugs,
    length(mapped_ids),
    drugs_with_moa,
    nrow(mechanisms),
    length(target_ids),
    human_targets,
    human_genes,
    gene_drug_cat
  )
)

fwrite(qc,
       file.path(DIR_PROCESSED, "ChEMBL_Preparation_QC.tsv"),
       sep="\t", quote=FALSE)

cat("\n============================================================\n")
cat("ChEMBL preparation completed\n")
cat("============================================================\n\n")
print(qc)

cat("\nOutputs:\n")
cat("  ", raw_dir, "/ChEMBL_drug_mapping.tsv\n", sep="")
cat("  ", raw_dir, "/ChEMBL_drug_mechanisms.tsv\n", sep="")
cat("  ", raw_dir, "/ChEMBL_targets.tsv\n", sep="")
cat("  ", raw_dir, "/ChEMBL_target_components.tsv\n", sep="")
cat("  ", raw_dir, "/ChEMBL_download_metadata.tsv\n", sep="")
cat("  ", file.path(DIR_PROCESSED, "ChEMBL_MoA_GeneDrugEvidence.tsv"), "\n", sep="")
cat("  ", file.path(DIR_PROCESSED, "ChEMBL_Preparation_QC.tsv"), "\n", sep="")

cat("\nNEXT STEP:\n")
cat("Inspect ChEMBL_Preparation_QC.tsv and unmapped drug rows before\n")
cat("running the direct mechanism-of-action benchmark.\n")
