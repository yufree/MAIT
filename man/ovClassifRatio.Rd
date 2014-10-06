\name{ovClassifRatioData}
\alias{ovClassifRatioData}
\title{Overall classification ratio extractor for MAIT objects}
\description{
Function \code{ovClassifRatioData} extracts the overall classification ratio for a \link{MAIT-class} object
}
\usage{

ovClassifRatioData(MAIT.object)

}

\arguments{
  \item{MAIT.object}{
A \link{MAIT-class} object
}
}
\value{
A list containing the overall classification ratio of the \link{MAIT-class} object for each classifier.
}
\seealso{
  \code{\link{MAIT-class}}  \code{\link{Validation}}}


\examples{
data(MAIT_sample)
MAIT<-spectralSigFeatures(MAIT,p.adj="fdr",parametric=TRUE)
MAIT <- Validation(Iterations = 20, trainSamples= 15, MAIT.object = MAIT)
ovClassifRatioData(MAIT)
}


\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}
