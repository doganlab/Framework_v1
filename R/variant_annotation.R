suppressPackageStartupMessages({
  library(data.table)
})

assert_columns <- function(dt, required, object_name = "data") {
  missing <- setdiff(required, names(dt))
  if (length(missing) > 0L) {
    stop(sprintf("%s is missing required columns: %s",
                 object_name, paste(missing, collapse = ", ")),
         call. = FALSE)
  }
  invisible(TRUE)
}

normalise_chr <- function(x) {
  x <- as.character(x)
  sub("^chr", "", x, ignore.case = TRUE)
}

prepare_vep_inputs <- function(instrument_file,
                               output_dir,
                               assembly = "GRCh37") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  instruments <- fread(instrument_file, na.strings = c("", "NA", "NaN"))
  required <- c("Metabolite", "SNP", "CHR", "BP",
                "EffectAllele", "OtherAllele")
  assert_columns(instruments, required, "InstrumentEvidence")

  instruments[, SNP := trimws(as.character(SNP))]
  instruments[, CHR := normalise_chr(CHR)]
  instruments[, BP := as.integer(BP)]
  instruments[, EffectAllele := toupper(trimws(as.character(EffectAllele)))]
  instruments[, OtherAllele := toupper(trimws(as.character(OtherAllele)))]

  invalid <- instruments[
    is.na(SNP) | SNP == "" |
      !grepl("^rs[0-9]+$", SNP) |
      is.na(CHR) | CHR == "" |
      is.na(BP) | BP < 1L
  ]

  valid <- instruments[
    !(is.na(SNP) | SNP == "" |
        !grepl("^rs[0-9]+$", SNP) |
        is.na(CHR) | CHR == "" |
        is.na(BP) | BP < 1L)
  ]

  # VEP ID input is deliberately used because EffectAllele/OtherAllele do not
  # establish REF/ALT orientation. VEP resolves rsIDs against the selected
  # assembly/cache, avoiding an incorrectly constructed VCF.
  vep_ids <- unique(valid[, .(SNP)])
  setorder(vep_ids, SNP)

  # Preserve every SNP-metabolite relation separately; one SNP can instrument
  # several metabolites and must not be collapsed at this stage.
  manifest_cols <- intersect(
    c("SNP", "Metabolite", "StudyAccession", "CHR", "BP",
      "EffectAllele", "OtherAllele", "EffectAlleleFrequency",
      "BetaExposure", "SEExposure", "PExposure", "FStatistic",
      "InstrumentStrengthClass", "WeakInstrument", "HarmonizationStatus",
      "Palindromic", "Ambiguous", "MHC", "SteigerDirection", "SourceFile"),
    names(valid)
  )
  manifest <- unique(valid[, ..manifest_cols])
  setorder(manifest, Metabolite, SNP)

  coordinate_check <- unique(valid[, .(
    SNP, CHR, BP, EffectAllele, OtherAllele
  )])
  coordinate_conflicts <- coordinate_check[, .(
    N_CoordinateAlleleCombinations = .N,
    CoordinateAlleleCombinations = paste(
      unique(paste(CHR, BP, EffectAllele, OtherAllele, sep = ":")),
      collapse = ";"
    )
  ), by = SNP][N_CoordinateAlleleCombinations > 1L]

  fwrite(vep_ids,
         file.path(output_dir, sprintf("vep_rsids_%s.txt", assembly)),
         col.names = FALSE, quote = FALSE)
  fwrite(manifest,
         file.path(output_dir, "VariantMetaboliteManifest.tsv"),
         sep = "\t", quote = FALSE, na = "NA")
  fwrite(invalid,
         file.path(output_dir, "InvalidVariantRecords.tsv"),
         sep = "\t", quote = FALSE, na = "NA")
  fwrite(coordinate_conflicts,
         file.path(output_dir, "VariantCoordinateConflicts.tsv"),
         sep = "\t", quote = FALSE, na = "NA")

  summary <- data.table(
    Metric = c(
      "Instrument rows",
      "Valid instrument rows",
      "Invalid instrument rows",
      "Unique rsIDs",
      "Unique metabolites",
      "rsIDs linked to multiple metabolites",
      "rsIDs with coordinate/allele conflicts",
      "Assembly requested"
    ),
    Value = c(
      nrow(instruments),
      nrow(valid),
      nrow(invalid),
      uniqueN(valid$SNP),
      uniqueN(valid$Metabolite),
      valid[, uniqueN(Metabolite), by = SNP][V1 > 1L, .N],
      nrow(coordinate_conflicts),
      assembly
    )
  )
  fwrite(summary,
         file.path(output_dir, "VEPInput_summary.tsv"),
         sep = "\t", quote = FALSE, na = "NA")

  invisible(list(
    instruments = valid,
    vep_ids = vep_ids,
    manifest = manifest,
    invalid = invalid,
    coordinate_conflicts = coordinate_conflicts,
    summary = summary
  ))
}
