\name{getScoresTable}
\alias{getScoresTable}
\title{Returns a list with the peak scores, masses, retention time and other information}
\description{
Function getScoresTable takes an \code{\link{MAIT-class}} object and returns a list with the scores of the features in the samples. Additionally, it returns the spectral ID of the peak and (optionally) a table containing the peak information (mass, retention time and annotation).
}
\usage{

getScoresTable(MAIT.object = NULL,
             getSpectra = TRUE,
             getExtendedTable = FALSE)
}
\arguments{
  \item{MAIT.object}{
A \link{MAIT-class} object where significant features have already been found.
}
\item{getSpectra}{
If it is set to TRUE, an element of the returned list will contain the spectra ID of each feature.
}

\item{getExtendedTable}{
If it is set to TRUE, an element of the returned list will contain a table with peak information (mass, retention time, annotation, intensity per sample)
}


}
\value{
 A list containing:
\itemize{
\item{
scores: The intensity of each feature per sample
}
\item{
spectraID: A numeric with the correspondence between peaks and spectral ID
}
\item{
extendedTable: a data frame containg detailed peak information (mass, retention time, annotation, intensity per sample).
}
}
}
\seealso{
\code{\link{spectralTStudent}} \code{\link{spectralAnova}}
}

\examples{
data(MAIT_sample)
MAIT<-spectralSigFeatures(MAIT,p.adj="fdr",parametric=TRUE)
getScoresTable(MAIT,getExtendedTable=TRUE)
}


\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}

