\name{t.tests.treatment.sign}
\alias{t.tests.treatment.sign}
\docType{data}
\title{
List of differentially expressed genes 
}
\description{
Consist of list of differentially expressed genes (DEG) with fold-change information i.e. up and down regulated genes between case and control.
}
\usage{data("t.tests.treatment.sign")}
\format{
  A data frame with 1820 observations on the following 8 variables.
  \describe{
    \item{\code{gene}}{a character vector}
    \item{\code{p.value}}{a numeric vector}
    \item{\code{p.value.fdr}}{a numeric vector}
    \item{\code{fc}}{a numeric vector}
    \item{\code{mean.A}}{a numeric vector}
    \item{\code{mean.B}}{a numeric vector}
    \item{\code{sd.A}}{a numeric vector}
    \item{\code{sd.B}}{a numeric vector}
  }
}
\examples{
data(t.tests.treatment.sign)
## maybe str(t.tests.treatment.sign) ;
}
\keyword{datasets}
