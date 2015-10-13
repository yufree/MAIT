\name{negAdducts}
\alias{negAdducts}
\title{Negative adducts table}
\description{
This table contains the adducts to be looked up for when the LC/MS polarisation mode was set to Negative. The layout of the table is that of the \code{CAMERA} adduct table.
}
\value{
A table having the fields:
\itemize{
\item{
ID: An ID number
}
\item{
name: The adduct name
}
\item{
nmol: Number of fragments in the adduct
}
\item{
charge: Electric charge of the adduct
}
\item{
massdiff: Mass difference in the fragment caused by the adduct
}
\item{
oidscore: Numeric relating the related clusters of ions
}
\item{
quasi: Binary value showing the validness of the annotation group
}
\item{
ips: Four values are possible (0.25,0.5,0.75,1) depending on the likelihood of the rule
}

}
}
\seealso{
  \code{\link{peakAnnotation}}}

\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}
\keyword{file}
