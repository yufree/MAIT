\name{LSDResults}
\alias{LSDResults}
\title{Extractor of the Fisher's LSD tests from a MAIT object}
\description{
Function \code{LSDResults} extracts the results of the LSD tests of a MAIT object as a matrix. 
}
\usage{

LSDResults(MAIT.object)

}

\arguments{
  \item{MAIT.object}{
A \link{MAIT-class} object
}
}
\value{
A matrix containing the results of the Fisher's LSD tests. For each row, equal letters mean that the groups are found to be equal in the test.
}
\seealso{
  \code{\link{MAIT-class}} \code{\link{FisherLSD}} \code{\link{spectralSigFeatures}} }

\examples{
data(MAIT_sample)
MAIT<-spectralSigFeatures(MAIT,p.adj="fdr",parametric=TRUE)
LSDResults(MAIT)
}

\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}
