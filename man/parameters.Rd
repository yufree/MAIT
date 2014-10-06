\name{parameters}
\alias{parameters}
\title{Extractor of the parameters used in the whole run from a MAIT object}
\description{
Function \code{parameters} extracts the slot link{MAIT.Parameters-class} of a \code{\link{MAIT-class}} object. This class contains all the parameters that have been used in the previous functions. Typing a summary of this object, a matrix version of the parameters is obtained.
}
\usage{

parameters(MAIT.object)

}

\arguments{
  \item{MAIT.object}{
A \link{MAIT-class} object
}
}
\value{
An object of the class MAIT.Parameters.
}
\seealso{
  \code{\link{MAIT-class}}}

\examples{
data(MAIT_sample)
MAIT<-spectralSigFeatures(MAIT,p.adj="fdr",parametric=TRUE)
MAIT <- Validation(Iterations = 20, trainSamples= 15, MAIT.object = MAIT)
parameters(MAIT)
}

\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}
