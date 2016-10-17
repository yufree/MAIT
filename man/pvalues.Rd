\name{pvalues}
\alias{pvalues}
\title{Pvalues extractor from a MAIT object}
\description{
Function \code{pvalues} extracts the pvalues contained in a \code{\link{MAIT-class}} object.
}
\usage{

pvalues(MAIT.object)

}

\arguments{
  \item{MAIT.object}{
A \link{MAIT-class} object
}
}
\value{
A numeric vector containing the pvalues of a \link{MAIT-class} object.
}
\seealso{
  \code{\link{MAIT-class}}  \code{\link{spectralSigFeatures}} }


\examples{
data(MAIT_sample)
MAIT<-spectralSigFeatures(MAIT,p.adj="fdr",parametric=TRUE)
pvalues(MAIT)
}

\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}
