\name{MAIT.FeatureData-class}
\Rdversion{1.1}
\docType{class}
\alias{MAIT.FeatureData-class}

\title{Class \code{"MAIT.FeatureData"}}
\description{
MAIT.FeatureData objects are used in the MAIT package to save the feature data.
}
\section{Slots}{
  \describe{
    \item{\code{scores}:}{Here it is saved the dataset obtained after applying the \link{peakAggregation} function.}
    \item{\code{featureID}:}{The ID numbers of all features are saved here.}
    \item{\code{featureSigID}:}{The ID numbers of the significant features are saved here.}
    \item{\code{LSDResults}:}{The results of performing a Fisher LSD test on each significant variable are saved in this slot.}
    \item{\code{models}:}{The model for each feature used to obtain the scores are saved in this slot.}
    \item{\code{pvalues}:}{In this slot are saved the pvalues of the features.}
    \item{\code{pvaluesCorrection}:}{The pvalues corrected by multiple test correction are saved here.}
    \item{\code{pcaModel}:}{PCA model generated using the function \code{\link{plotPCA}}.}
    \item{\code{plsModel}:}{PCA model generated using the function \code{\link{plotPLS}}.}
    \item{\code{masses}:}{Masses used as an input for the function \code{\link{MAITbuilder}}.}
    \item{\code{rt}:}{Retention time values used as an input for the function \code{\link{MAITbuilder}}.}
    \item{\code{extendedTable}:}{Dataframe containing the information regarding masses, retention time values, intensity and spectra IDs
          passes as an input for the function \code{\link{MAITbuilder}}}


  }
}
\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}
\seealso{
  \code{\link{MAIT-class}}       \code{\link{xsAnnotate}}}
\keyword{classes}


                  