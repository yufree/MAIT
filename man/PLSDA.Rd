\name{PLSDA}
\alias{PLSDA}
\title{Applies PLSDA to the provided data}
\description{
Function \code{PLSDA} performs Fisher's LSD tests on the data using the package plsgenomics
}
\usage{

PLSDA(Xtrain, 
      Ytrain, 
      Xtest = NULL, 
      ncomp,
      nruncv = 0, 
      alpha = 2/3, 
      priors = NULL) 

}

\arguments{
  \item{Xtrain}{
A numerical matrix containing the data
}

  \item{Ytrain}{
A factor vector containing the class labels of the samples
}

  \item{Xtest}{
A numerical matrix containing the data whose class is to be predicted.
}

  \item{ncomp}{
Number of components to build the PCA model
}

  \item{nruncv}{
Number of cross-validation iterations to be performed for the choice of the number of latent components
}

  \item{alpha}{
The proportion of Obesrvations to be included in the training set at each cross-validation iteration
}

  \item{priors}{
The class priors to be used for linear discriminant analysis. If unspecified, the class proportions in the training set are used.
}
}

\value{
A list containing the output of function pls.regression, the predicted class for the Xtest dataset and the number of components used.
}

\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}
\keyword{internal}