\name{spectralTStudent}
\alias{spectralTStudent}
\title{Extract significant features from a MAIT object for two classes}
\description{
Function spectralTStudent takes a \code{\link{MAIT-class}} object and obtains which of the variables are significant given a p-value threshold when there only are two classes in the raw data. The parameters of the significant features can be printed to an output table (TRUE by default).
}
\usage{

spectralTStudent(MAIT.object = NULL, 
		 pvalue = 0.05, 
		 p.adj = "none", 
                 printCSVfile = TRUE)
}
\arguments{
  \item{MAIT.object}{
A \link{MAIT-class} object where function \link{peakAggregation} has already been applied. The output of the function is going to be an update of the same \link{MAIT-class} object.
}
  \item{pvalue}{
P-value threshold. Variables having a p-value lower than this value is considered as a significant variable.
}

  \item{p.adj}{
Character with the name of the posthoc method to be applied to correct the pvalues. The supported methods are that of the \link{p.adjust} function
}
  \item{printCSVfile}{
Set to TRUE if an output table has to be produced. See function \link{sigPeaksTable} for more information.
}

}
\value{
  A \link{MAIT-class} object containing the significant features of the scores slot of \link{MAIT-class} object used as an input.
}
\seealso{
 \code{\link{spectralSigFeatures}} \code{\link{sigPeaksTable}}
}
\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}
\keyword{internal}
