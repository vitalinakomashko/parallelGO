% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/not_exported.R
\name{return_error_result}
\alias{return_error_result}
\title{Output a data frame for sets that generated an error or a warning}
\usage{
return_error_result(error_warning, input_genes, input_label)
}
\arguments{
\item{error_warning}{A condition object.}

\item{input_genes}{Character vector with the list of genes for which GO
enrichment is run.}

\item{input_label}{Character string with a label for the list of genes.}
}
\value{
Data frame with 1 row.
}
\description{
\code{return_error_result} creates a one-row data frame with the same columns
as the regular output of get_ontology. The column \strong{count} is
represented by the number of genes, the column \strong{term} provides the
message, and the column \strong{set_label} provides the label of the gene set
for which a condition was encountered. This is an internal function.
}
\seealso{
\code{\link{run_go}}, \code{\link{remove_errors}},
\code{\link{get_ontology}}.
}
