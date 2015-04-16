\name{pvaluesCorrection}
\alias{pvaluesCorrection}
\title{P-values correction extractor from a MAIT object}
\description{
Function \code{pvaluesCorrection} returns a character showing wheter some multiple testing correction has been perfomed on the p-values.
}
\usage{

pvaluesCorrection(MAIT.object)

}

\arguments{
  \item{MAIT.object}{
A \link{MAIT-class} object
}
}
\value{
The output is a character whose values could be "None" if no p-value correction has been performed or "Bonferroni" if Bonferroni multiple test correction was selected when function \link{spectralSigFeatures} was applied.
}
\seealso{
  \code{\link{MAIT-class}}  \code{\link{spectralSigFeatures}} }

\examples{
data(MAIT_sample)
MAIT<-spectralSigFeatures(MAIT,p.adj="fdr",parametric=TRUE)
pvaluesCorrection(MAIT)
}


\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}
