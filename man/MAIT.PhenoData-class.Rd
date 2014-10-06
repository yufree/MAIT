\name{MAIT.PhenoData-class}
\Rdversion{1.1}
\docType{class}
\alias{MAIT.PhenoData-class}

\title{Class \code{"MAIT.PhenoData"}}
\description{
MAIT.PhenoData objects are used in the MAIT package to save the phenotype data.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("MAIT.PhenoData", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
	The information related to the classes present in the data is stored in this slot. It has three different extra slots:
  \describe{
    \item{\code{classes}:}{It contains the name of the classes in the data. It can be quickly accessed by using the function \code{\link{classes}}}
    \item{\code{classNum}:}{Vector showing the number of samples belonging to each class. It can be quickly accessed by using the function \code{\link{classNum}}}
    \item{\code{resultsPath}:}{In this slot is saved the direction where the project is saved. This means that all the output tables and files of the MAIT object are going to be stored in that directory. It can be quickly accessed by using the function \code{\link{resultsPath}}}
  }
}
\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}
\seealso{
   \code{\link{MAIT-class}}      \code{\link{xsAnnotate}}}
\keyword{classes}
