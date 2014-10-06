\name{selectPLScomp}
\alias{selectPLScomp}
\title{Looks for the optimum number of components to be considered for the PLSDA}
\description{
Function \code{selectPLScomp} finds the optimum number of components to be used by the Partial Least Squares and linear Discriminant Algorithm (PLSDA).
}
\usage{

selectPLScomp(data,
         class,
         max.comp)

}

\arguments{
  \item{data}{
A numerical matrix containing the data
}

  \item{class}{
Vector containing the class label of each sample.
}

  \item{max.comp}{
Maximum number of components to be considered.
}


}

\value{
A numeric value of the optimal number of components to be considered in a PLSDA run.
}
\seealso{
  \code{\link{Validation}}}

\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}
\keyword{internal}