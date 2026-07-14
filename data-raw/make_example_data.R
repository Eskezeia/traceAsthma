# Generates the small, fully synthetic example dataset bundled with the
# package (inst/extdata/*.csv and data/example_cohort.rda), used in the
# vignette and in package examples. Five patients, two SNPs, two CpGs,
# four target genes, three TFs -- deliberately tiny so the full pipeline
# runs in seconds and every intermediate number can be hand-checked, as
# used throughout the design discussion for this package.

set.seed(42)

patient_ids <- sprintf("P%03d", 1:5)

# --- genotype dosage (2 SNPs) ---
dosage <- matrix(
  c(0, 1, 2, 1, 0,
    2, 1, 0, 1, 2),
  nrow = 5, dimnames = list(patient_ids, c("rs_TSLP", "rs_IL33"))
)

# --- methylation (2 CpGs, beta-like values in [0,1]) ---
methylation <- matrix(
  c(0.72, 0.35, 0.58, 0.66, 0.44,
    0.61, 0.42, 0.39, 0.55, 0.48),
  nrow = 2, byrow = TRUE,
  dimnames = list(c("cg_IL13", "cg_GATA3"), patient_ids)
)

# --- expression (4 target genes, TPM-like) ---
expression <- matrix(
  c(250,  80, 190, 150, 110,
    180,  70, 120,  90,  85,
    310,  90, 260, 210, 140,
    160,  60, 130, 100,  95),
  nrow = 4, byrow = TRUE,
  dimnames = list(c("IL4", "IL5", "IL13", "GATA3"), patient_ids)
)

# --- TF-target regulon (3 TFs) ---
regulon_db <- data.frame(
  tf = c("GATA3", "GATA3", "GATA3", "STAT6", "STAT6", "NFKB1"),
  target = c("IL4", "IL5", "IL13", "IL13", "IL4", "IL4"),
  mor = c(1, 1, 1, 1, 1, -1),
  confidence = c("A", "A", "A", "B", "B", "C"),
  stringsAsFactors = FALSE
)

# --- clinical / outcome ---
clinical <- data.frame(
  patient_id = patient_ids,
  asthma = c(1, 0, 1, 1, 0),
  age = c(34, 41, 28, 52, 37),
  sex = c("F", "F", "M", "F", "M"),
  bmi = c(24.1, 27.8, 22.6, 31.2, 25.9),
  smoking = c("never", "former", "never", "current", "never"),
  stringsAsFactors = FALSE
)

covariates <- data.frame(
  age = clinical$age,
  sex = ifelse(clinical$sex == "F", 1, 0),
  pc1 = c(0.02, -0.01, 0.03, -0.02, 0.01),
  pc2 = c(-0.01, 0.02, 0.00, 0.01, -0.02),
  row.names = patient_ids
)

eqtm_pairs <- data.frame(
  cpg = c("cg_IL13", "cg_IL13", "cg_GATA3", "cg_GATA3"),
  gene = c("IL13", "GATA3", "IL4", "IL5"),
  distance_bp = c(500, NA, NA, NA),
  stringsAsFactors = FALSE
)
same_chr <- c(TRUE, FALSE, FALSE, FALSE)

example_cohort <- list(
  dosage = dosage,
  methylation = methylation,
  expression = expression,
  regulon_db = regulon_db,
  clinical = clinical,
  covariates = covariates,
  eqtm_pairs = eqtm_pairs,
  same_chr = same_chr
)

usethis_available <- requireNamespace("usethis", quietly = TRUE)
if (usethis_available) {
  usethis::use_data(example_cohort, overwrite = TRUE)
} else {
  if (!dir.exists("data")) dir.create("data")
  save(example_cohort, file = "data/example_cohort.rda")
}

# Also write flat CSVs to inst/extdata for the Shiny app demo / non-R users
write.csv(methylation, "inst/extdata/example_methylation.csv")
write.csv(expression, "inst/extdata/example_expression.csv")
write.csv(clinical, "inst/extdata/example_clinical.csv", row.names = FALSE)
