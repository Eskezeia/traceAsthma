\name{import_data}
\alias{import_data}
\alias{export_data}
\title{Import and export omics/clinical data in multiple real-world file formats}
\usage{
import_data(path, format = NULL, as_matrix = TRUE, sheet = 1, con = NULL, ...)

export_data(x, path, format = NULL, ...)
}
\arguments{
\item{path}{File path (CSV, TSV, \code{.xlsx}/\code{.xls}, \code{.rds}, \code{.parquet}) or,
for \code{format = "database"}, a table name to query via \code{con}.}
\item{format}{One of \code{"csv"}, \code{"tsv"}, \code{"excel"}, \code{"rds"}, \code{"parquet"},
\code{"database"}. If \code{NULL} (default), inferred from \code{path}'s extension.}
\item{as_matrix}{Logical; if \code{TRUE} (default), the first column is used
as row names and the remainder coerced to a numeric matrix. Set \code{FALSE}
for clinical/covariate tables, which should remain data frames.}
\item{sheet}{For \code{format = "excel"}, the sheet name or index to read. Default 1.}
\item{con}{For \code{format = "database"}, an open \code{DBI} connection.}
\item{...}{Additional arguments passed to the underlying reader/writer.}
\item{x}{A data frame or matrix, for \code{export_data()}.}
}
\value{
\code{import_data()} returns a numeric matrix (if \code{as_matrix = TRUE}) or a
data frame. \code{export_data()} invisibly returns \code{path}.
}
\description{
A single entry point for reading and writing real-world data in the
shapes \code{traceAsthma}'s pipeline functions expect, so that transitioning
from the bundled simulated/example data to a real dataset requires
changing only the data-loading call, not any downstream analysis code.
Supports CSV, TSV, Excel, RDS, Parquet, and (for import) database
connections via \code{DBI}.
}
\examples{
sim <- simulate_trace_asthma_cohort(n_subjects = 10, n_cpgs = 5,
                                     n_genes = 5, n_tfs = 2, seed = 1)
tmp <- tempfile(fileext = ".csv")
export_data(sim$methylation, tmp)
meth <- import_data(tmp)
dim(meth)
}
