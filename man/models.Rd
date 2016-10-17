\name{models}
\alias{models}
\title{Model extractor from a MAIT object}
\description{
Function \code{models} extracts the models of a MAIT object as a list. 
}
\usage{

models(MAIT.object)

}

\arguments{
  \item{MAIT.object}{
A \link{MAIT-class} object
}
}
\value{
A list containing the models of the MAIT object.
}
\seealso{
  \code{\link{MAIT-class}}}

\examples{
data(MAIT_sample)
MAIT<-spectralSigFeatures(MAIT,p.adj="fdr",parametric=TRUE)
models(MAIT)
}

\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}
