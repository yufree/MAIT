\name{annotateBiotransf}
\alias{annotateBiotransf}
\title{Single Biotransformation Annotator}
\description{
Function \code{annotateBiotransf} annotates a spectrum with an already detected Biotransformation.
}
\usage{

annotateBiotransf(biotransf, diffIndex, spectrum, sigPeaksTable,biotransformationsTable)

}

\arguments{
	  \item{biotransf}{
The already detected biotransformation using the \code{\link{inBetween}} function
}
  \item{diffIndex}{
A numeric pointer to the other peak involved in the biotransformation inside the spectrum.
}
  \item{spectrum}{
The spectrum to be annotated.
}
  \item{sigPeaksTable}{
A dataframe obtained from running the \code{\link{sigPeaksTable}} function.
}
  \item{biotransformationsTable}{
Table of biotransformations either read from the bioTable argument or the default MAIT table.
}

}
\value{
A vector containg the masses of the peaks involved in the biotransformation, their retention time, their annotation and their indices in the signPeaksTable.
}
\seealso{
  \code{\link{Biotransformations}}}

\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}

\keyword{internal}