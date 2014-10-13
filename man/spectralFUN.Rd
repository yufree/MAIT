\name{spectralFUN}
\alias{spectralFUN}
\title{Extract significant features from a MAIT object using a user-defined test}
\description{
Function spectralFUN takes an \code{\link{MAIT-class}} object and obtains which of the variables are significant given a p-value threshold following a user-defined statistical test. The parameters of the significant features can ve printed to an output table (TRUE by default).
}
\usage{

spectralFUN(pvalue=0.05,
            p.adj="none",
            MAIT.object=NULL,
            printCSVfile=TRUE,
            test.fun=NULL,
            namefun=NULL)}

\arguments{

  \item{pvalue}{
P-value threshold. Variables having a p-value lower than this value is considered as a significant variable.
}

  \item{p.adj}{
Post-hoc method to be used to correct the p-values.
}

  \item{MAIT.object}{
A \link{MAIT-class} object where function \link{peakAggregation} has already been applied. The output of the function is going to be an update of the same \link{MAIT-class} object.
}

  \item{printCSVfile}{
Set to TRUE if an output table has to be produced. See function \link{sigPeaksTable} for more information.
}

  \item{test.fun}{
Function containing the statistical test to be applied on each feature. The function should be designed to correct just one feature as the function will aplly this correction to all the features in the MAIT.object.
}

  \item{namefun}{
Character with the name of the test. This name will appear in the MAITparameters table and in the summary of the MAIT object.
}


}
\value{
  A \link{MAIT-class} object containing the significant features of the scores slot of \link{MAIT-class} object used as an input.
}
\seealso{
  \code{\link{MAIT-class}} \code{\link{peakAggregation}} \code{\link{sigPeaksTable}}
}
\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}
\keyword{internal}