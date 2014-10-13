\name{MAIT-class}
\Rdversion{1.1}
\docType{class}
\alias{MAIT-class}
\alias{model,MAIT-method}
\alias{summary,MAIT-method}
\alias{scores,MAIT-method}
\alias{loadings,MAIT-method}

\title{Class \code{"MAIT"}}
\description{
MAIT class objects are used in the MAIT package to perform the analysis and statistical calculus of LC/MS data. It has 5 main slots: FeatureInfo, RawData, Validation, PhenoData and FeatureData
}
\section{Slots}{
  \describe{
    \item{\code{FeatureInfo}:}{Object of class \code{\link{MAIT.FeatureInfo-class}}}
    \item{\code{RawData}:}{Object of class \code{\link{MAIT.RawData-class}}}
    \item{\code{Validation}:}{Object of class \code{\link{MAIT.Validation-class}}}
    \item{\code{PhenoData}:}{Object of class \code{\link{MAIT.PhenoData-class}} }
    \item{\code{FeatureData}:}{Object of class \code{\link{MAIT.FeatureData-class}} }
  }
}
\section{Methods}{
  \describe{
    \item{summary}{\code{signature(object = "MAIT")}: This function show a summary of the workflow results performed so far including the classification results and the parameters used. }
    \item{model}{\code{signature(object = "MAIT")}: returns the model for either the PCA, PLS models when functions \code{\link{plotPCA}} or \code{\link{plotPLS}} have been already respectively launched. }
    \item{scores}{\code{signature(object = "MAIT")}: Retrieves the scores from a MAIT object}
    \item{loadings}{\code{signature(object = "MAIT")}: Retrieves the loadings from a MAIT object}

	 }
}
\seealso{
  \code{\link{xsAnnotate}}}

\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}
\keyword{classes}
