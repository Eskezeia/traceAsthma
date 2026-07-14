\name{data_dictionary}
\alias{data_dictionary}
\alias{create_data_templates}
\title{Data dictionary and blank templates for mapping real datasets}
\usage{
data_dictionary(object = NULL)

create_data_templates(dir, include_example_row = TRUE)
}
\arguments{
\item{object}{Optional character, one of \code{"methylation"}, \code{"expression"},
\code{"dosage"}, \code{"clinical"}, \code{"regulon_db"}, \code{"eqtm_pairs"}; if supplied,
returns only that object's entries.}
\item{dir}{Directory to write templates into. Created if it doesn't exist.}
\item{include_example_row}{Logical; if \code{TRUE} (default), each template
includes one illustrative row of correctly formatted example values.}
}
\value{
\code{data_dictionary()} returns a data frame with columns \code{object}, \code{field},
\code{type}, \code{units}, \code{coding_scheme}, \code{required}, \code{description}.
\code{create_data_templates()} invisibly returns a character vector of the
file paths written.
}
\description{
\code{data_dictionary()} documents every data object and field
\code{traceAsthma}'s pipeline functions expect -- the reference for mapping a
real-world dataset onto package requirements, independent of whatever
variable names or coding schemes an upstream source system uses.
\code{create_data_templates()} writes blank (or example-row) CSV templates
matching the dictionary, so a user or collaborator preparing a real
dataset has a concrete file to fill in and validate against.
}
\examples{
dd <- data_dictionary("clinical")
print(dd, row.names = FALSE)

tmp <- tempfile()
create_data_templates(tmp)
}
