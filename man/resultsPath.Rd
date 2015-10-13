\name{resultsPath}
\alias{resultsPath}
\title{Retrieves the folder where the results are saved for a MAIT object}
\description{
Function \code{resultsPath} returns the folder where the plots and tables are saved for a \code{\link{MAIT-class}} object
}  
\usage{

resultsPath(MAIT.object)

}

\arguments{
  \item{MAIT.object}{
A \link{MAIT-class} object
}
}
\value{
A character showing where the plots and tables have been stored.
}
\seealso{
  \code{\link{MAIT-class}}}

\examples{
data(MAIT_sample)
MAIT<-spectralSigFeatures(MAIT,p.adj="fdr",parametric=TRUE)
resultsPath(MAIT)
}

\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}
