#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
source(file.path(dirname(script_path), "..", "R", "bootstrap.R"))
source(file.path(PROJECT_ROOT, "R", "io.R"))
source(file.path(PROJECT_ROOT, "R", "utils.R"))
source(file.path(PROJECT_ROOT, "R", "pharmacometabolomics.R"))

check_file_exists(RAW_COLOC_FILE)
check_file_exists(NATMED_WORKBOOK_FILE)

coloc <- data.table::fread(RAW_COLOC_FILE)
metabolites <- read_natmed_metabolite_dictionary(NATMED_WORKBOOK_FILE, NATMED_METABOLITE_SHEET)
drugs <- read_natmed_drug_categories(NATMED_WORKBOOK_FILE, NATMED_DRUG_CATEGORY_SHEET)
assoc <- read_natmed_m1_m2_associations(NATMED_WORKBOOK_FILE, NATMED_ASSOCIATION_SHEET)
mapping <- map_coloc_to_natmed(coloc, metabolites)
matched <- mapping[MatchStatus == "matched"]
assoc_matched <- assoc[NatMedAbbreviation %in% matched$NatMedAbbreviation]

observed <- data.table(
  Metric = c(
    "CAD_coloc_traits",
    "NatMed_traits",
    "Matched_traits",
    "NatMed_drug_categories",
    "Suppl_T4_rows",
    "Matched_T4_rows",
    "Matched_M2_significant_rows",
    "Matched_M2_significant_drug_categories"
  ),
  Observed = c(
    uniqueN(coloc$metabolite),
    uniqueN(metabolites$NatMedAbbreviation),
    uniqueN(matched$metabolite),
    uniqueN(drugs$DrugCategory),
    nrow(assoc),
    nrow(assoc_matched),
    assoc_matched[M2Significant %in% TRUE, .N],
    uniqueN(assoc_matched[M2Significant %in% TRUE, DrugCategory])
  ),
  Expected_for_current_files = c(
    177L,
    150L,
    102L,
    87L,
    13050L,
    8874L,
    1260L,
    42L
  )
)

observed[, PASS := Observed == Expected_for_current_files]
print(observed)

# A few specific structural checks that matter for the framework.
extra <- data.table(
  Test = c(
    "XXL_VLDL_P maps to XXL.VLDL.P",
    "All matched traits have one NatMed abbreviation",
    "All Suppl_T4 drug categories map to Suppl_T1",
    "Clopidogrel has M2-significant associations in matched traits"
  ),
  PASS = c(
    mapping[metabolite == "XXL_VLDL_P", unique(NatMedAbbreviation)][1] == "XXL.VLDL.P",
    mapping[MatchStatus == "matched", .(N = uniqueN(NatMedAbbreviation)), by = metabolite][, all(N == 1)],
    all(assoc$DrugCategory %in% drugs$DrugCategory),
    assoc_matched[DrugCategory == "B_Clopidogrel" & M2Significant %in% TRUE, .N] > 0
  )
)
print(extra)

if (any(!observed$PASS) || any(!extra$PASS)) {
  stop("Nature Medicine smoke test FAILED.")
}

cat("\nNature Medicine smoke test PASSED.\n")
cat("Note: the current source files yield 102, not 103, matched Nightingale traits.\n")
cat("Acetone is not in the NatMed 150-trait panel; Acetate (Ace) is a different trait.\n\n")
