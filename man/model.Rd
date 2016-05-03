\name{model}
\alias{model}
\title{Model extractor for either PCA or PLS models}
\description{
Function \code{model} returns the model for either the PCA, PLS models when functions \code{\link{plotPCA}} or \code{\link{plotPLS}} have been already respectively launched. 
}
\usage{

model(x,type)

}

\arguments{
  \item{x}{
A \link{MAIT-class} object
}
  \item{type}{
A character whose value should be "PCA" or "PLS" depending on which loading vectors are wanted.
}

}
\value{
The PCA or PLS model
}
\seealso{
  \code{\link{plotPCA}} or \code{\link{plotPLS}}}

\examples{
data(MAIT_sample)
MAIT<-spectralSigFeatures(MAIT,p.adj="fdr",parametric=TRUE)
MAIT<-plotPCA(MAIT,plot3d=FALSE)
model(MAIT,type="PCA")

MAIT<-plotPLS(MAIT,plot3d=FALSE)
model(MAIT,type="PLS")
}


\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}
