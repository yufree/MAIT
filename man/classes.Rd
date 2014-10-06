\name{classes}
\alias{classes}
\title{Class names extractor from a MAIT object}
\description{
Function \code{classes} extracts the class names of a link{MAIT-class} object as a vector. 
}
\usage{

classes(MAIT.object)

}

\arguments{
  \item{MAIT.object}{
A \link{MAIT-class} object
}
}

\value{
A character vector containing the class names of the \link{MAIT-class} object
}
\seealso{
  \code{\link{MAIT-class}}}

\examples{
data(MAIT_sample)
MAIT
classes(MAIT)
}

\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}
