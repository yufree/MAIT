\name{MAIT.FeatureInfo-class}
\Rdversion{1.1}
\docType{class}
\alias{MAIT.FeatureInfo-class}

\title{Class \code{"MAIT.FeatureInfo"}}
\description{
MAIT.FeatureInfo objects are used in the MAIT package to save the data related to the information of the features.
}
\section{Slots}{
	It stores information related to the features. It contains three extra slots:
  \describe{
    \item{\code{biotransformations}:}{Biotransformations found when function \link{Biotransformations} is launched. }
    \item{\code{peakAgMethod}:}{In this slot is stored the table created by the function \link{identifyMetabolites}. It can be retrieved quickly in R by using the function \link{metaboliteTable} }
    \item{\code{metaboliteTable}:}{Peak Aggregation Method used when function \link{peakAggregation} is launched.}
  }
}
\section{Methods}{
No methods defined with class "MAIT.FeatureInfo" in the signature.
}
\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
    \code{\link{MAIT-class}}     \code{\link{xsAnnotate}}}
\keyword{classes}
