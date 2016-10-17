\name{sigPeaksTable}
\alias{sigPeaksTable}
\title{Build a table of the information related to the significant features contained in a MAIT object}
\description{
Function sigPeaksTable takes an \code{\link{MAIT-class}} object containing significant feature information and builds a table with the information related to these features.
}
\usage{

sigPeaksTable(MAIT.object=NULL,
             printCSVfile=FALSE,
             extendedTable = TRUE,
             printAnnotation=TRUE)
}
\arguments{
  \item{MAIT.object}{
A \link{MAIT-class} object where significant features have already been found.
}
\item{printCSVfile}{
Set to TRUE if an output table has to be produced. The table should be found in (working directory)/(project directory)Tables/significativeFeatures.csv. 
}

\item{extendedTable}{
Set to TRUE the table created by the peak external data is used.
}

\item{printAnnotation}{
Set to TRUE The peak annotation is provided in the output table
}


}
\value{
 A table containing:
\itemize{
\item{
First column (mz): Peak mass
}
\item{
Second column(mzmin): Minimum peak mass of the peak group.
}
\item{
Third column(mzmax): Maximum peak mass of the peak group.
}
\item{
Fourth column(rt): Peak retention time (in minutes).
}
\item{
Fifth column(rtmin): Minimum peak retention time of the peak group.
}
\item{
Sixth column(rtmax): Maximum peak retention time of the peak group.
}
\item{ 
Seventh column(npeaks): Number of samples where the peak has been detected.
}
\item{
The columns from the nineth to the column labeled "isotopes" contain number of class samples where the peak has been detected and the intensities of the peak among samples. 
}
\item{
The isotopes column shows if the peak has been identified as a possible isotope.
}
\item{
The adduct column shows which kind of adduct could the peak be.
}
\item{
The column labeled pcgroup contains the spectral ID of the peak.
}
\item{
The P.adjust column contains the corrected peak p-value using post-hoc methods.
}
\item{
The p column shows the peak p-value with no multiple test correction.
}
\item{
The Fisher column shows the Fisher test results for the peak. Each of the letters separated by the character "_" corresponds to a class value. Classes having the same letters are indistinguible whereas those having different letters are statistically different clases.
}
\item{
The last columns contain the mean and median values for each feature
}
}
}
\seealso{
\code{\link{spectralTStudent}} \code{\link{spectralAnova}}
}

\examples{
data(MAIT_sample)
MAIT<-spectralSigFeatures(MAIT,p.adj="fdr",parametric=TRUE)
head(sigPeaksTable(MAIT))
}


\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}
\keyword{file}
