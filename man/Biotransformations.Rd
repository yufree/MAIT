\name{Biotransformations}
\alias{Biotransformations}
\title{Biotransformations of the significant data contained in a MAIT object are identified.}
\description{
This function takes a \link{MAIT-class} object having significant features already detected and looks up for biotransformations between them. \pkg{MAIT} has a default biotransformation table that will be used if no other table is specified via the bioTable input parameter.
}
\usage{

Biotransformations(MAIT.object = NULL, 
                   peakPrecision = 0.005,
                   bioTable = NULL, 
                   adductTable = NULL,
                   adductAnnotation = FALSE)
}
\arguments{
  \item{MAIT.object}{
A \link{MAIT-class} object where significant features have already been found.
}
  \item{peakPrecision}{
Maximum difference between the peak masses differences and the values shown in bioTable to be considered as a biotransformation. As default the value is 0.005 Da.
}
  \item{bioTable}{
Table containing the biotransformations to be looked for in the signData input. By default it is taken the \pkg{MAIT-class} biotransformations table.
}

  \item{adductTable}{
Table containing the adducts to be looked for in the signData input. By default it is taken the \pkg{MAIT-class} positive adducts table. If this argument is set to "negAdducts", then the default table for negative adducts is taken instead. It is possible to use a user-defined adduct table
}
  \item{adductAnnotation}{
If it is set to TRUE, both adduct and Biotransformations annotation stages are performed.
}


}
\value{
A \link{MAIT-class} object with the updated biotransformations slot
}
\seealso{
  \code{\link{spectralAnova}}     \code{\link{spectralTStudent}}}

\examples{
data(MAIT_sample)
MAIT<-spectralSigFeatures(MAIT,p.adj="fdr",parametric=TRUE)
MAIT<-Biotransformations(MAIT.object = MAIT, peakPrecision = 0.005)
MAIT@FeatureInfo@biotransformations; #Detected Biotransformations
}



\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}
\keyword{file}
