\name{format_nested_model_table}
\alias{format_nested_model_table}
\alias{format_coefficient_table}
\alias{format_mediation_table}
\alias{export_table}
\title{Publication-standard (Nature-style) table export}
\usage{
format_nested_model_table(nested_comparison, auc_ci = NULL)

format_coefficient_table(risk_model, variable_labels = NULL)

format_mediation_table(sem_result)

export_table(table, path)
}
\arguments{
\item{nested_comparison}{Output of \code{compare_nested_models()}.}
\item{auc_ci}{Optional data frame with columns \code{model}, \code{ci_low},
\code{ci_high}.}
\item{risk_model}{A fitted \code{glm} object (e.g. from \code{fit_risk_model()}).}
\item{variable_labels}{Optional named character vector mapping raw model
term names to publication-ready labels, e.g.
\code{c(eqtm = "Methylation eQTM score (per SD)")}. Strongly recommended
so the exported table never shows a bare, non-descriptive variable name.}
\item{sem_result}{Output of \code{fit_mediation_sem()}.}
\item{table}{A \code{traceAsthma_table} object (from a \code{format_*_table()}
function).}
\item{path}{Output file path (\code{.csv} or \code{.xlsx}).}
}
\value{
The \code{format_*_table()} functions return a \code{traceAsthma_table}
object (a data frame with a \code{footnote} attribute); print it directly
or pass it to \code{export_table()}. \code{export_table()} invisibly
returns the output path.
}
\description{
Produces tables matching Nature-family journal conventions: fully
spelled-out column headers, predictors and models named directly (never
"M1"-"M5" alone), consistent decimal formatting, and a footnote
explaining every abbreviation and statistical convention used --
required content for a Nature-style table that a reader could not get
from the column headers alone. \code{export_table()} writes the table as
CSV with a companion \verb{_footnote.txt} file (dependency-free default),
or as a single-sheet \code{.xlsx} with the footnote appended below the
table if the \code{writexl} package is installed.
}
\examples{
\dontrun{
tbl <- format_nested_model_table(nested_comparison)
print(tbl)
export_table(tbl, "table1_nested_models.csv")
}
}
