\name{rawData}
\alias{rawData}
\title{Raw data extractor from a MAIT object}
\description{
Function \code{rawData} extracts the raw data used to build the \code{\link{MAIT-class}} object
}
\usage{

rawData(MAIT.object)

}

\arguments{
  \item{MAIT.object}{
A \link{MAIT-class} object
}
}
\value{
A list containing either a xcmsSet or a \link{xsAnnotate} object.
}
\seealso{
   \code{\link{xsAnnotate-class}}}

\examples{
data(MAIT_sample)
MAIT
rawData(MAIT)
}

\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}
