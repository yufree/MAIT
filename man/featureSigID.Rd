\name{featureSigID}
\alias{featureSigID}
\title{Feature statistically significant ID extractor from a MAIT object}
\description{
Function \code{featureSigID} extracts the vector index of the feature IDs of a MAIT object that have been found significant through function \link{spectralSigFeatures}.
}
\usage{

featureSigID(MAIT.object)

}

\arguments{
  \item{MAIT.object}{
A \link{MAIT-class} object
}
}
\value{
A numeric vector containing the statistically significant feature IDs of the MAIT object.
}
\seealso{
  \code{\link{MAIT-class}} \code{\link{spectralSigFeatures}}}

\examples{
data(MAIT_sample)
MAIT<-spectralSigFeatures(MAIT,p.adj="fdr",parametric=TRUE)
featureSigID(MAIT)
featureID(MAIT)[featureSigID(MAIT)] #Significant spectra IDs
}


\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}
