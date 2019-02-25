\name{human_symbol}
\alias{human_symbol}
\docType{data}
\title{
 List of human gene symbols grouped by set label.
}
\description{
Small test dataset with human gene symbols created for demonstrating
package functionality.
}
\usage{data("human_symbol")}
\format{
  A data frame with 280 observations on the following 2 variables.
  \describe{
    \item{\code{id}}{gene symbols}
    \item{\code{set_label}}{gene set labels}
  }
}
\source{
internal
}
\examples{
data(human_symbol)
str(human_symbol)
}
\keyword{datasets}
