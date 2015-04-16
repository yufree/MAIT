\name{classifRatioClasses}
\alias{classifRatioClasses}
\title{Class classification ratio extractor from a MAIT object}
\description{
Function \code{classifRatioClasses} extracts the class classification ratio of a MAIT object as a matrix. 
}
\usage{

classifRatioClasses(MAIT.object)

}

\arguments{
  \item{MAIT.object}{
A \link{MAIT-class} object where function \link{Validation} has already been launched successfully.
}
}
\value{
A matrix containing the classification ratio for each class, classifier and iteration.
}
\seealso{
  \code{\link{MAIT-class}}}

\examples{
data(MAIT_sample)
MAIT<-spectralSigFeatures(MAIT,p.adj="fdr",parametric=TRUE)
MAIT <- Validation(Iterations = 20, trainSamples= 15, MAIT.object = MAIT)
classifRatioClasses(MAIT)
}


\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}
