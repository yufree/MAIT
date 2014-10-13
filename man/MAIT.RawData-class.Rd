\name{MAIT.RawData-class}
\Rdversion{1.1}
\docType{class}
\alias{MAIT.RawData-class}

\title{Class \code{"MAIT.RawData"}}
\description{
MAIT.RawData objects are used in the MAIT package to save the data related to the information of the features.
}
\section{Slots}{
	This class contains information of the raw data and the parameters of the whole analysis. It has two slots:
  \describe{
    \item{\code{parameters}:}{All the parameters of the analysis are saved in this slot. It can be obtained as a matrix in R by typing summary(parameters(MAIT.object))}
    \item{\code{data}:}{This slot contains either the xcmsSet-class or the \link{xsAnnotate} object, depending if the function \link{peakAnnotation} has already been launched}
  }
}
\section{Methods}{
No methods defined with class "MAIT.RawData" in the signature.
}
\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}

\seealso{
  \code{\link{MAIT-class}}     \code{\link{xsAnnotate}}}

\keyword{classes}
