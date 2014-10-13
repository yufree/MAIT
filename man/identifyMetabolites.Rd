\name{identifyMetabolites}
\alias{identifyMetabolites}
\title{Metabolite identifier}
\description{
Takes a MAIT object and performs the metabolite search for the significant features
}
\usage{

identifyMetabolites(MAIT.object=NULL,
		    peakTolerance=0.005,
                    database=NULL,
                    polarity="positive",
                    printCSVfile=TRUE)
}
\arguments{
  \item{MAIT.object}{
A \link{MAIT-class} object where significant features have already been found.
}
  \item{peakTolerance}{
Maximum difference between the peak masses differences and the values shown in the database to be considered as a match. As default the value is 0.005 Da.
}
  \item{database}{
User-defined input table. If it is set to NULL, the default \pkg{MAIT} database is selected to perform the metabolite identification.
}
  \item{polarity}{
Character parameter that can be set to "positive" or "negative" depending on the polarity in which the samples were taken.
}
\item{printCSVfile}{
Set to TRUE if an output table has to be produced. The table should be found in (working directory)/Tables/SearchTable.csv. 
}


}
\value{
An output table is stored in the folder (working directory)/Tables/SearchTable.csv if \code{printCSVfile} is set to TRUE. More info at \link{metaboliteTable}

}
\seealso{
  \code{\link{Biotransformations}}  \code{\link{spectralSigFeatures}}    }

\examples{
data(MAIT_sample)
MAIT<-spectralSigFeatures(MAIT,p.adj="fdr",parametric=TRUE)
MAIT<-Biotransformations(MAIT.object = MAIT, peakPrecision = 0.005)
MAIT <- identifyMetabolites(MAIT.object = MAIT, peakTolerance = 0.005,polarity="positive")
}


\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}
\keyword{file}
