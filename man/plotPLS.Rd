\name{plotPLS}
\alias{plotPLS}
\title{2D and 3D PLS scoreplots from a MAIT object}
\description{
This function takes a MAIT-class object containing information of the significant features in the data and performs 2D scoreplots (PC1 vs PC2, PC2 vs PC3 and PC1 vs PC3) saved as png files. Additionally it also performs an interactive 3D PLS scoreplot.
}
\usage{

plotPLS (MAIT.object=NULL,
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
If the number of components in the PLS is found to be three or more, three different PLS scoreplots are printed in three png files. One using PC1 vs PC2, another with PC1 vs PC3 and the last one with PC2 vs PC3. If the number of components is less than three, all the possible plots of these three are created. The files will be stored in the directory (working directory)/PLS_Scoreplots. Moreover, an interactive 3D PLS scoreplot is also generated through function \link{PCAplot3d}.
}
\seealso{
  \code{\link{spectralAnova}}     \code{\link{spectralTStudent}}}

\examples{
data(MAIT_sample)
MAIT<-spectralSigFeatures(MAIT,p.adj="fdr",parametric=TRUE)
MAIT<-plotPLS(MAIT,plot3d=FALSE)
}


\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}
\keyword{file}
