\name{method}
\alias{method}
\title{Peak Aggregation Method Used}
\description{
Function \code{method} returns the name of the peak aggregation method used on a \link{MAIT-class} object.
}
\usage{

method(object)

}

\arguments{
  \item{object}{
A \link{MAIT-class} object
}

}
\value{
A character with the peak aggregation method
}
\seealso{
  \code{\link{peakAggregation}}}

\examples{
data(MAIT_sample)
method(MAIT)
}

\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}
\keyword{internal}