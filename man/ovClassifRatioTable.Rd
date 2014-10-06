\name{ovClassifRatioTable}
\alias{ovClassifRatioTable}
\title{Overall classification table extractor for MAIT objects}
\description{
Function \code{ovClassifRatioData} extracts the overall classification table for a \code{\link{MAIT-class}} object
}
\usage{

ovClassifRatioTable(MAIT.object)

}

\arguments{
  \item{MAIT.object}{
A \link{MAIT-class} object
}
}
\value{
A list containing the overall classification table of the MAIT object for each classifier showing the mean value and their standard error.
}
\seealso{
  \code{\link{MAIT-class}}  \code{\link{Validation}}}

\examples{
data(MAIT_sample)
MAIT<-spectralSigFeatures(MAIT,p.adj="fdr",parametric=TRUE)
MAIT <- Validation(Iterations = 20, trainSamples= 15, MAIT.object = MAIT)
ovClassifRatioTable(MAIT)
}


\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}
