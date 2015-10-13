\name{plotHeatmap}
\alias{plotHeatmap}
\title{Builds ten heatmaps with different p-values and clustering distances}
\description{
This function takes a MAIT object containing information of the significant features in the data and plots 10 heatmaps. 5 diferent p-values (0.05, 0.01, 0.001, 1e-4 and 1e-5) and two clustering distances (euclidean and pearson) are used.
}
\usage{

plotHeatmap(MAIT.object=NULL)
}
\arguments{
  \item{MAIT.object}{
A \link{MAIT-class} object where significant features have already been found.
}
}
\value{
10 different heatmaps using 5 p-values (0.05, 0.01, 0.001, 1e-4 and 1e-5) and two clustering distances (euclidean and pearson) are created. The plots will be stored as png files in a folder called (working directory)/Heatmaps
}
\seealso{
  \code{\link{spectralAnova}}     \code{\link{spectralTStudent}}}

\examples{
data(MAIT_sample)
MAIT<-spectralSigFeatures(MAIT,p.adj="fdr",parametric=TRUE)
MAIT<-plotHeatmap(MAIT)
}



\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}
\keyword{function}
