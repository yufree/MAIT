\name{MAIT.Validation-class}
\Rdversion{1.1}
\docType{class}
\alias{MAIT.Validation-class}

\title{Class \code{"MAIT.Validation"}}
\description{
MAIT.Validation objects are used in the MAIT package to save the validation results obtained from the classification run.
}
\section{Slots}{
	The information related to the run of the function \link{Validation} is saved here. It contains three lots:

  \describe{
    \item{\code{ovClassifRatioTable}:}{Summary table showing the overall classification ratios for each of the three classifiers. It can be quickly gathered by using the function \code{ovClassifRatioTable}}
    \item{\code{ovClassifRatioData}:}{All the data corresponding to the overall classification ratios. It can be quickly gathered by using the function \code{\link{ovClassifRatioData}}}
    \item{\code{classifRatioClasses}:}{All the data corresponding to the classification ratios per class. It can be quickly gathered by using the function \code{\link{classifRatioClasses}}}
  }
}
\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}

\seealso{
   \code{\link{MAIT-class}}      \code{\link{xsAnnotate}}}
\keyword{classes}
