#' Data dictionary: expected variables, types, units, and coding schemes
#'
#' Returns a structured table documenting every data object and field
#' `traceAsthma`'s pipeline functions expect -- the single reference for
#' mapping a real-world dataset onto package requirements, independent of
#' whatever variable names/coding schemes an upstream source system uses.
#' Print it, filter it by `object`, or use it to drive a custom mapping
#' script from your source data to the shapes this package expects.
#'
#' @param object Optional character, one of `"methylation"`,
#'   `"expression"`, `"dosage"`, `"clinical"`, `"regulon_db"`,
#'   `"eqtm_pairs"`; if supplied, returns only that object's entries.
#' @return A data frame with columns `object`, `field`, `type`, `units`,
#'   `coding_scheme`, `required`, `description`.
#' @export
#' @examples
#' dd <- data_dictionary()
#' print(dd[dd$object == "clinical", ], row.names = FALSE)
data_dictionary <- function(object = NULL) {
  dd <- data.frame(
    object = c(
      rep("methylation", 3), rep("expression", 3), rep("dosage", 3),
      rep("clinical", 9), rep("regulon_db", 4), rep("eqtm_pairs", 3)
    ),
    field = c(
      # methylation
      "row names", "column names", "matrix values",
      # expression
      "row names", "column names", "matrix values",
      # dosage
      "row names", "column names", "matrix values",
      # clinical
      "patient_id", "asthma", "age", "sex", "bmi", "smoking",
      "ancestry_pc1", "ancestry_pc2", "(additional ancestry_pc* columns)",
      # regulon_db
      "tf", "target", "mor", "confidence",
      # eqtm_pairs
      "cpg", "gene", "distance_bp"
    ),
    type = c(
      "character (CpG probe ID)", "character (subject ID)", "numeric",
      "character (gene ID/symbol)", "character (subject ID)", "numeric",
      "character (SNP/variant ID)", "character (subject ID)", "numeric",
      "character", "integer", "numeric", "numeric (0/1)", "numeric",
      "factor/character", "numeric", "numeric", "numeric",
      "character (TF gene symbol)", "character (target gene ID/symbol)",
      "numeric (+1/-1)", "character (A/B/C)",
      "character (CpG probe ID)", "character (gene ID/symbol)", "numeric or NA"
    ),
    units = c(
      "-", "-", "beta value, proportion methylated in [0,1] (or M-value, see allow_mvalues)",
      "-", "-", "raw or normalized counts (TPM/CPM/VST), non-negative",
      "-", "-", "allele dosage, 0-2 copies of the effect allele",
      "-", "-", "years", "-", "kg/m^2", "-", "-", "-", "-",
      "-", "-", "-", "-",
      "-", "-", "base pairs, distance from CpG to gene TSS"
    ),
    coding_scheme = c(
      "e.g. \"cg00000029\" (Illumina probe ID)", "must match clinical$patient_id / other matrices' column names", "-",
      "e.g. gene symbol \"IL13\" or Ensembl ID", "must match clinical$patient_id / other matrices' column names", "-",
      "e.g. \"rs123456\"", "must match clinical$patient_id / other matrices' column names", "-",
      "unique subject identifier, must match all matrix column names", "0 = control, 1 = case (asthma)",
      "-", "0 = male, 1 = female (see validate_clinical() for auto-detection of F/M text codes)",
      "-", "levels: \"never\", \"former\", \"current\"",
      "genetic ancestry principal component 1 (from PCA on genome-wide genotypes)",
      "genetic ancestry principal component 2", "additional PCs as needed, named ancestry_pc3, ancestry_pc4, ...",
      "e.g. \"GATA3\"", "e.g. \"IL4\"", "+1 = activation, -1 = repression",
      "DoRothEA-style confidence tier, A = highest",
      "must match a row name in the methylation matrix", "must match a row name in the expression matrix",
      "NA or > cis window (e.g. 1e6) for trans pairs on the same chromosome; see same_chr for different-chromosome pairs"
    ),
    required = c(
      TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
      TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
      TRUE, TRUE, TRUE, FALSE,
      TRUE, TRUE, FALSE
    ),
    description = c(
      "One row per CpG probe.", "One column per subject; order must match all other objects (see validate_cohort_alignment()).",
      "Methylation beta value for that CpG in that subject.",
      "One row per gene.", "One column per subject; order must match all other objects.",
      "Expression level for that gene in that subject.",
      "One row per genetic variant.", "One column per subject; order must match all other objects.",
      "Imputed or hard-called dosage of the effect allele.",
      "Primary key linking the clinical table to matrix column names.",
      "Outcome variable for all risk models.", "Age at data collection/enrollment.",
      "Biological sex, used as a covariate.", "Body mass index, used as a covariate if available.",
      "Smoking status, used as a covariate if available.",
      "Used to adjust for population stratification; omit if single-ancestry cohort.",
      "As above.", "traceAsthma auto-detects any column matching \"^ancestry_pc\".",
      "Transcription factor.", "Target gene of that TF.",
      "Direction of regulation for that TF-target pair.",
      "Evidence tier; used to filter regulons via build_regulons(confidence_classes = ...).",
      "Candidate CpG for eQTM discovery.", "Candidate target gene for eQTM discovery.",
      "Used to classify the pair as cis or trans in discover_eqtm()."
    ),
    stringsAsFactors = FALSE
  )
  if (!is.null(object)) dd <- dd[dd$object == object, ]
  dd
}

#' Write blank CSV templates for every input object, for real-data mapping
#'
#' Generates empty (header-only, or header + one illustrative example row)
#' CSV templates matching [data_dictionary()], so a user or collaborator
#' preparing a real dataset has a concrete file to fill in and validate
#' against, rather than reverse-engineering the expected format from
#' function documentation alone.
#'
#' @param dir Directory to write templates into. Created if it doesn't exist.
#' @param include_example_row Logical; if `TRUE` (default), each template
#'   includes one illustrative row of correctly-formatted example values
#'   (clearly not real data) in addition to the header.
#' @return Invisibly, a character vector of the file paths written.
#' @export
#' @examples
#' tmp <- tempfile()
#' paths <- create_data_templates(tmp)
#' basename(paths)
create_data_templates <- function(dir, include_example_row = TRUE) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)

  clinical_template <- data.frame(
    patient_id = if (include_example_row) "P0001" else character(0),
    asthma = if (include_example_row) 1 else integer(0),
    age = if (include_example_row) 34 else numeric(0),
    sex = if (include_example_row) 0 else numeric(0),
    bmi = if (include_example_row) 24.1 else numeric(0),
    smoking = if (include_example_row) "never" else character(0),
    ancestry_pc1 = if (include_example_row) 0.012 else numeric(0),
    ancestry_pc2 = if (include_example_row) -0.004 else numeric(0)
  )

  regulon_template <- data.frame(
    tf = if (include_example_row) "GATA3" else character(0),
    target = if (include_example_row) "IL4" else character(0),
    mor = if (include_example_row) 1 else numeric(0),
    confidence = if (include_example_row) "A" else character(0)
  )

  eqtm_pairs_template <- data.frame(
    cpg = if (include_example_row) "cg00000029" else character(0),
    gene = if (include_example_row) "IL13" else character(0),
    distance_bp = if (include_example_row) 50000 else numeric(0)
  )

  methylation_template <- if (include_example_row) {
    m <- matrix(c(0.42, 0.87), nrow = 2, dimnames = list(c("cg00000029", "cg00000108"), "P0001"))
    data.frame(cpg_id = rownames(m), m, check.names = FALSE)
  } else {
    data.frame(cpg_id = character(0))
  }

  paths <- c(
    clinical = file.path(dir, "clinical_template.csv"),
    regulon_db = file.path(dir, "regulon_db_template.csv"),
    eqtm_pairs = file.path(dir, "eqtm_pairs_template.csv"),
    methylation = file.path(dir, "methylation_template.csv")
  )
  utils::write.csv(clinical_template, paths["clinical"], row.names = FALSE)
  utils::write.csv(regulon_template, paths["regulon_db"], row.names = FALSE)
  utils::write.csv(eqtm_pairs_template, paths["eqtm_pairs"], row.names = FALSE)
  utils::write.csv(methylation_template, paths["methylation"], row.names = FALSE)

  message("Templates written to ", dir, ":\n  ", paste(basename(paths), collapse = "\n  "),
          "\n(Expression and dosage matrices follow the same layout as the methylation template: ",
          "first column = feature ID, remaining columns = one per subject.)")
  invisible(paths)
}
