\name{selectK}
\alias{selectK}
\title{Looks for the optimum number of nearest neighbours to be considered for the KNN}
\description{
Function \code{selectK} finds the optimum number of nearest neighbours for the K-Nearest Neighbours (KNN) algorithm.
}
\usage{

selectK(data,
         class,
         max.k)

}

\arguments{
  \item{data}{
A numerical matrix containing the data
}

  \item{class}{
Vector containing the class label of each sample.
}

  \item{max.k}{
Maximum number of nearest neighbours to be considered.
}


}

\value{
A numeric value of the optimal number of neighbours to be considered in a KNN run.
}
\seealso{
  \code{\link{Validation}}}

\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}
\keyword{internal}