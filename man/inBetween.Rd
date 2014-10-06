\name{inBetween}
\alias{inBetween}
\title{Checks if a peak mass value is in a certain mass allowance window.}
\description{
Function \code{inBetween} extracts the mass peaks of a certain spectrum provided a dataframe where the spectrum labels are in a column called pcgroup.
}
\usage{

inBetween(testValue,biotRange)

}

\arguments{
  \item{testValue}{
The peak mass value to be checked
}
  \item{biotRange}{
A matrix containing two numerical columns and each row refers to a certain neutral mass loss. The first column should contain the lower value (neutral mass value minus the peak allowance window) and the second column should have the higher value (neutral mass value plus the peak allowance window)
}
}
\value{
The rows of the biotRange table where possible neutral losses have been detected.
}
\seealso{
  \code{\link{Biotransformations}}}

\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}
\keyword{internal}