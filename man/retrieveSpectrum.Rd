\name{retrieveSpectrum}
\alias{retrieveSpectrum}
\title{Extractor of the mass peaks corresponding to a certain spectrum}
\description{
Function \code{retrieveSpectrum} extracts the mass peaks of a certain spectrum provided a dataframe where the spectrum labels are in a column called pcgroup.
}
\usage{

retrieveSpectrum(spectrumNumber,sigPeaksTable)

}

\arguments{
  \item{spectrumNumber}{
The spectrum ID number whose peaks we want to retrieve.
}
  \item{sigPeaksTable}{
A dataframe containing the peak data in rows. There should be a column called pcgroup containing the spectra correspondence for all the peaks and the firts column should contain the peak masses.
}
}
\value{
A numeric vector containing the peak masses of the queried spectrum.
}
\seealso{
  \code{\link{Biotransformations}}}

\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}
\keyword{internal}