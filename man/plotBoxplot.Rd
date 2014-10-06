\name{plotBoxplot}
\alias{plotBoxplot}
\title{Prints a png file for each of the significant peak/spectra present in the input}
\description{
This function takes a \link{MAIT-class} object containing information related to the significant features and plots a boxplot for each significant feature (peak or spectra).
}
\usage{

plotBoxplot(MAIT.object=NULL)
}
\arguments{
  \item{MAIT.object}{
A \link{MAIT-class} object where significant features have already been found.
}

}
\value{
A boxplot is stored as a png file for each of the significant features (peak or spectra). The files will be stored in the directory (working directory)/Boxplots
}
\seealso{
  \code{\link{spectralAnova}} \code{\link{spectralTStudent}}}

\examples{
data(MAIT_sample)
MAIT<-spectralSigFeatures(MAIT,p.adj="fdr",parametric=TRUE)
MAIT<-plotBoxplot(MAIT)
}




\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}
\keyword{function}
