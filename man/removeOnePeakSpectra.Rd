\name{removeOnePeakSpectra}
\alias{removeOnePeakSpectra}
\title{Removes those spectra having just one peak}
\description{
Function \code{removeOnePeakSpectra} removes the spectra having just one peak
}
\usage{

removeOnePeakSpectra(data, 
      idGroup) 

}

\arguments{
  \item{data}{
A numerical matrix containing the peak data
}

  \item{idGroup}{
A numeric vector containing the spectra id number of the peaks
}
}


\value{
A peak data set without the one-peak spectra.
}
\seealso{
  \code{\link{peakAggregation}}}

\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}
\keyword{internal}