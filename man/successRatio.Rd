\name{successRatio}
\alias{successRatio}
\title{Extracts the success ratio of a truth table}
\description{
Function \code{successRatio} extracts the success ratio (weighted ratio of samples correctly classified vs total samples) for each class and overall. The value is weighted to take into account the possible different sample number between classes.
}
\usage{

successRatio(classes,
             tt,
             ClassWeights)

}

\arguments{
  \item{tt}{
The truth table from which the success ratio should have to be extracted.
}

  \item{classes}{
Vector containing the class label of each sample.
}

  \item{ClassWeights}{
Vector containing the weights of each class
}


}

\value{
A numeric value showing the ratio of the samples that are correctly associated according to their real classes.
}
\seealso{
  \code{\link{Validation}}}

\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}
\keyword{internal}