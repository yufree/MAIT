\name{PCAplot3d}
\alias{PCAplot3d}
\title{3D PCA scoreplots}
\description{
This function takes three sets of coordinates and builds a 3D scoreplot using the package \pkg{rgl}
}
\usage{

PCAplot3d (z,
           x,
           y,
           cols,
           axes=TRUE,
           new=TRUE)
}
\arguments{
  \item{z}{
A numerical vector containing the values for the z-axis.
}
  \item{y}{
A numerical vector containing the values for the y-axis.
}
  \item{x}{
A numerical vector containing the values for the x-axis.
}

  \item{axes}{
Boolean parameter. Set to TRUE if axes should be plotted.
}
  \item{new}{
Boolean parameter. Set to TRUE if a new rgl plot should be created.
}
  \item{cols}{
Character vector containing the colors for each sample.
}



}
\value{
A 3D interactive plot is created using the package \pkg{rgl}
}
\seealso{
  \code{\link{plotPCA}}}

\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}
\keyword{internal}
