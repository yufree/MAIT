\name{classNum}
\alias{classNum}
\title{Sample number extractor for each class from a MAIT object}
\description{
Function \code{classNum} extracts the number of samples belonging to each class of a MAIT object as a vector. 
}
\usage{

classNum(MAIT.object)

}

\arguments{
  \item{MAIT.object}{
A \link{MAIT-class} object
}
}
\value{
A numeric vector containing the number of samples for each class of the MAIT object. The order of the classes correspond to that of the output of function \link{classes}.
}
\seealso{
  \code{\link{MAIT-class}}}

\examples{
data(MAIT_sample)
MAIT
classNum(MAIT)
}


\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}