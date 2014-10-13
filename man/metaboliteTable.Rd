\name{metaboliteTable}
\alias{metaboliteTable}
\title{Metabolite table generator}
\description{
Takes a \code{\link{MAIT-class}} object and builds a table with the information related to the significant features and their possible identifications.
}
\usage{

metaboliteTable(MAIT.object,
                printCSVfile = FALSE)}
\arguments{
  \item{MAIT.object}{
A \link{MAIT-class} object where significant features have already been found.
}
  \item{printCSVfile}{
A boolean parameter. Set to TRUE if a csv file should be written with the metabolite table.
}}


\value{
An output table is stored in the folder (working directory)/Tables/SearchTable.csv having the fields:
\itemize{
\item{
First column: search ID number.
}
\item{
Second column (mz): Peak mass.
}
\item{
Third column(rt): Peak retention time (in minutes).
}
\item{
The columns from the third to the column labeled "p.adj" contain number of class samples where the peak has been detected and the intensities of the peak among samples. 
}
\item{
The P.adjust column contains the corrected peak p-value using bonferroni.
}
\item{
The p column shows the peak p-value with no multiple test correction.
}
\item{
The Fisher column shows the FIsher test results for the peak. Each of the letters separated by the character "_" corresponds to a class value. Classes having the same letters are indistinguible whereas those having different letters are statistically different clases.
}
\item{
The isotopes column shows if the peak has been identified as a possible isotope.
}
\item{
The adduct column shows which kind of adduct or biotransformation could the peak be.
}
\item{
Column Name contains the name of the possible metabolite identification for the peak.
}
\item{
The column labeled spectra contains the spectral ID of the peak.
}
\item{
Column Biofluid shows if the identified search is stored as a biofluid in the input database or not.
}
\item{
The column ENTRY shows the database name of the entry for the metabolite.
}
}
}
\seealso{
  \code{\link{identifyMetabolites}}  \code{\link{spectralAnova}}     \code{\link{spectralTStudent}}}

\examples{
data(MAIT_sample)
MAIT<-spectralSigFeatures(MAIT,p.adj="fdr",parametric=TRUE)
MAIT<-Biotransformations(MAIT.object = MAIT, peakPrecision = 0.005)
MAIT <- identifyMetabolites(MAIT.object = MAIT, peakTolerance = 0.005,polarity="positive")
head(metaboliteTable(MAIT))
}


\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}
\keyword{file}
