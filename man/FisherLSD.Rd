\name{FisherLSD}
\alias{FisherLSD}
\title{Performs Fisher's LSD tests on the provided data}
\description{
Function \code{FisherLSD} performs Fisher's LSD tests on the data using the package agricolae.
}
\usage{

FisherLSD(data,
          classes,
          index,
          DFerror,
          MSerror,
          numClasses
          )

}

\arguments{
  \item{data}{
A numerical matrix containing the data
}

  \item{classes}{
A character vector containing the class names of the samples present in the data. This vector must have the same length as the number of samples present in the argument data.
}

  \item{index}{
Numerical value to choose a subset of the data on which the LSD tests is going to be performed.
}

  \item{DFerror}{
Degrees of freedom of the model
}

  \item{MSerror}{
Means square error of the model
}
  \item{numClasses}{
Numerical parameter corresponding to the number of classes present in the data.
}

}

\value{
A list containing the class names, the group where each class belongs according to the LSD test and the value of their means.
}


\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}
\keyword{internal}