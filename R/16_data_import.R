#' Import an omics or clinical data matrix/table from a file or connection
#'
#' A single entry point for reading real-world data into the shapes
#' `traceAsthma`'s pipeline functions expect, so that transitioning from
#' the bundled simulated/example data to a real dataset requires changing
#' only the data-loading call, not any downstream analysis code. Format is
#' auto-detected from the file extension unless `format` is supplied
#' explicitly.
#'
#' @param path File path (CSV, TSV, `.xlsx`/`.xls`, `.rds`, `.parquet`) or,
#'   for `format = "database"`, a table name to query via `con`.
#' @param format One of `"csv"`, `"tsv"`, `"excel"`, `"rds"`, `"parquet"`,
#'   `"database"`. If `NULL` (default), inferred from `path`'s extension.
#' @param as_matrix Logical; if `TRUE` (default for omics data types),
#'   the first column is used as row names and the remainder coerced to a
#'   numeric matrix (the shape [discover_eqtm()], [infer_tf_activity()],
#'   etc. expect for `methylation`/`expression`/`dosage`). Set `FALSE` for
#'   clinical/covariate tables, which should remain data frames.
#' @param sheet For `format = "excel"`, the sheet name or index to read
#'   (passed to `readxl::read_excel()`). Default 1 (first sheet).
#' @param con For `format = "database"`, an open `DBI` connection.
#' @param ... Additional arguments passed to the underlying reader
#'   (`utils::read.csv()`, `readxl::read_excel()`, `arrow::read_parquet()`,
#'   or `DBI::dbReadTable()`).
#' @return A numeric matrix (if `as_matrix = TRUE`) or a data frame.
#' @export
#' @examples
#' # Round-trip a simulated cohort through CSV to demonstrate the workflow
#' # that later applies unchanged to a real file.
#' sim <- simulate_trace_asthma_cohort(n_subjects = 20, n_cpgs = 15,
#'                                      n_genes = 10, n_tfs = 3, seed = 1)
#' tmp <- tempfile(fileext = ".csv")
#' utils::write.csv(sim$methylation, tmp)
#' meth <- import_data(tmp, format = "csv", as_matrix = TRUE)
#' dim(meth)
import_data <- function(path, format = NULL, as_matrix = TRUE, sheet = 1,
                         con = NULL, ...) {
  if (is.null(format)) {
    ext <- tolower(tools::file_ext(path))
    format <- switch(ext,
      csv = "csv", tsv = "tsv", txt = "tsv",
      xlsx = "excel", xls = "excel",
      rds = "rds", parquet = "parquet",
      stop("Could not infer format from extension '.", ext, "'. Supply `format` explicitly ",
           "(one of \"csv\", \"tsv\", \"excel\", \"rds\", \"parquet\", \"database\").", call. = FALSE)
    )
  }

  df <- switch(format,
    csv = utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, ...),
    tsv = utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, ...),
    excel = {
      requirePkg("readxl", "reading Excel files")
      as.data.frame(readxl::read_excel(path, sheet = sheet, ...))
    },
    rds = {
      obj <- readRDS(path)
      if (is.matrix(obj)) return(obj)
      as.data.frame(obj)
    },
    parquet = {
      requirePkg("arrow", "reading Parquet files")
      as.data.frame(arrow::read_parquet(path, ...))
    },
    database = {
      if (is.null(con)) stop("`con` (an open DBI connection) is required for format = \"database\".", call. = FALSE)
      requirePkg("DBI", "reading from a database connection")
      DBI::dbReadTable(con, path, ...)
    },
    stop("Unrecognized format '", format, "'.", call. = FALSE)
  )

  if (!as_matrix) return(df)

  rn <- as.character(df[[1]])
  mat <- as.matrix(df[, -1, drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- rn
  mat
}

#' Export a data frame or matrix to a file, matching [import_data()]'s formats
#'
#' @param x A data frame or matrix.
#' @param path Output file path; format inferred from extension unless
#'   `format` is supplied.
#' @param format One of `"csv"`, `"tsv"`, `"excel"`, `"rds"`, `"parquet"`.
#' @param ... Passed to the underlying writer.
#' @return Invisibly, `path`.
#' @export
export_data <- function(x, path, format = NULL, ...) {
  if (is.null(format)) {
    ext <- tolower(tools::file_ext(path))
    format <- switch(ext, csv = "csv", tsv = "tsv", txt = "tsv",
                      xlsx = "excel", rds = "rds", parquet = "parquet",
                      stop("Could not infer format from extension '.", ext, "'.", call. = FALSE))
  }
  if (is.matrix(x)) {
    x <- data.frame(id = rownames(x), x, check.names = FALSE, stringsAsFactors = FALSE)
  }
  switch(format,
    csv = utils::write.csv(x, path, row.names = FALSE, ...),
    tsv = utils::write.table(x, path, sep = "\t", row.names = FALSE, ...),
    excel = { requirePkg("writexl", "writing Excel files"); writexl::write_xlsx(x, path, ...) },
    rds = saveRDS(x, path, ...),
    parquet = { requirePkg("arrow", "writing Parquet files"); arrow::write_parquet(x, path, ...) },
    stop("Unrecognized format '", format, "'.", call. = FALSE)
  )
  invisible(path)
}
