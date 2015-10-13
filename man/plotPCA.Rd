\name{plotPCA}
\alias{plotPCA}
\title{2D and 3D PCA scoreplots from a MAIT object}
\description{
This function takes a MAIT-class object containing information of the significant features in the data and performs 2D scoreplots (PC1 vs PC2, PC2 vs PC3 and PC1 vs PC3) saved as png files. Additionally it also performs an interactive 3D PCA scoreplot.
}
\usage{

plotPCA (MAIT.object=NULL,
         Log=FALSE,
	 center=TRUE,
	 scale=TRUE,
	 plot3d=TRUE)
}
\arguments{
  \item{MAIT.object}{
A \link{MAIT-class} object where significant features have already been found.
}
  \item{Log}{
Set to TRUE if the data should be plotted using the logarithm of the intensity.
}
  \item{center}{
Set to TRUE if the data should be centered around its mean. See \code{\link{scale}}.
}

  \item{scale}{
Set to TRUE if the data should be scaled. See \code{\link{scale}}.
}

  \item{plot3d}{
Boolean set to TRUE if a 3D PCA scoreplot should be plot.
}


}
\value{
Three different PCA scoreplots are printed in three png files. One using PC1 vs PC2, another with PC1 vs PC3 and the last one with PC2 vs PC3. The files will be stored in the directory (working directory)/PCA_Scoreplots. Moreover, an interactive 3D PCA scoreplot is also generated through function \link{PCAplot3d}.
}
\seealso{
  \code{\link{spectralAnova}}     \code{\link{spectralTStudent}}}

\examples{
data(MAIT_sample)
MAIT<-spectralSigFeatures(MAIT,p.adj="fdr",parametric=TRUE)
MAIT<-plotPCA(MAIT,plot3d=FALSE)
}

\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}
\keyword{file}
