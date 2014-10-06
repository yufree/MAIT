\name{spectralWelch}
\alias{spectralWelch}
\title{Extract significant features from a MAIT object}
\description{
Function spectralWelch takes an \code{\link{MAIT-class}} object and obtains which of the variables are significant given a p-value threshold following a Welch test. The parameters of the significant features can ve printed to an output table (TRUE by default).
}
\usage{

spectralWelch(MAIT.object = NULL,
              pvalue = 0.05, 
              p.adj="none", 
              printCSVfile = TRUE)
}
\arguments{
  \item{MAIT.object}{
A \link{MAIT-class} object where function \link{peakAggregation} has already been applied. The output of the function is going to be an update of the same \link{MAIT-class} object.
}
  \item{pvalue}{
P-value threshold. Variables having a p-value lower than this value is considered as a significant variable.
}

  \item{p.adj}{
Post-hoc method to be used to correct the p-values.
}

  \item{printCSVfile}{
Set to TRUE if an output table has to be produced. See function \link{sigPeaksTable} for more information.
}

}
\value{
  A \link{MAIT-class} object containing the significant features of the scores slot of \link{MAIT-class} object used as an input.
}
\seealso{
  \code{\link{MAIT-class}} \code{\link{peakAggregation}} \code{\link{sigPeaksTable}}
}
\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}
\keyword{internal}
