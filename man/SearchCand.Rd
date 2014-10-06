\name{SearchCand}
\alias{SearchCand}
\title{Peak search function into a database}
\description{
Function \code{SearchCand} looks up for a peak into a database
}
\usage{

SearchCand(candidate,
           dataBase,
           peakTolerance)

}

\arguments{
  \item{candidate}{
The mass of the peak to be looked up into the database
}
  \item{dataBase}{
The table where the database to be used is saved. The function is build to use databases with the same layout as the MAIT's database. This database can be accessed by typing data(MAITtables) and Database.
}
  \item{peakTolerance}{
Maximum difference between the peak masses differences and the values shown in the database to be considered as a match.

}

}
\value{
A matrix containing all the possible hits for that peak candidate
}
\seealso{
  \code{\link{MAIT-class}}  \code{\link{identifyMetabolites}} }

\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}
\keyword{internal}