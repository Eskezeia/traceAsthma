\name{list_dependencies}
\alias{list_dependencies}
\title{List traceAsthma's dependencies and what each is used for}
\usage{
list_dependencies()
}
\value{
Invisibly, a data frame with columns \code{package}, \code{type}, \code{stage}.
}
\description{
Prints a table of every non-base package \code{traceAsthma} uses, whether
it is required or optional-per-stage (all are \code{Suggests}), and which
pipeline stage needs it, so a user can install only what their planned
analysis requires. QC, eQTM discovery, regulon filtering/enrichment, the
final logistic risk model, all publication table exporters
(\code{format_*_table}), and all figure functions (\code{plot_*}) work
with base R alone and need nothing from this list.
}
\examples{
list_dependencies()
}
