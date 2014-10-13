\name{peakAggregation}
\alias{peakAggregation}
\title{Performs a peak aggregation procedure to the rawData of a MAIT object}
\description{
  peakAggregation function applies a peak aggregation technique to the data of a \link{MAIT-class} object. Several aggregation techniques are available (see methods below).
}
\usage{

peakAggregation(MAIT.object=NULL,
                method="None",
                clases=NULL,
                samples=NULL,
                PCAscale=FALSE,
                PCAcenter=FALSE,
                scale=FALSE,
                signVariables=NULL,
                RemoveOnePeakSpectra=FALSE,
                printCSVfile=TRUE)
}
\arguments{
  \item{MAIT.object}{
A \link{MAIT-class} object where function \link{peakAnnotation} has already been applied. The output of the function is going to be an update of the same \link{MAIT-class} object.
}
  \item{method}{
Chosen method to perform the dimensionality reduction using the non-free \pkg{pagR} package:
- If it is set to "None", no reduction is performed and the spectral peaks are taken as variables. This is the default method. 
- If it is set to "Mean", the intensity mean value over each sample is taken and used as spectral intensity.
- If it is set to "PCA", the first scores vector of a principal components analysis (PCA) decomposition is used as spectral intensity.
- If it is set to "NMF", the first scores vector of a non-negative matrix factorization (NMF) is used as spectral intensity.
- If it is set to "Single", the spectral peak having the highest intensity mean value over samples among all the spectral peaks is used as spectral intensity.
}
  \item{clases}{
Parameter to explicitly define the clases of the future spectralData object. If it is set to NULL this value is taken from the annotatedPeaks input. 

}
  \item{samples}{
If the spectralData object has to include just a subset of the annotatedPeaks' samples, this input must be the vector having the wanted sample's IDs.
 }
  \item{PCAscale}{
If method="PCA" and PCAscale is set to TRUE, then the data is scaled following the \code{\link{prcomp}} function. If it is set to TRUE, scale input is ignored.
  }
\item{PCAcenter}{
If method="PCA" and PCAscale is set to TRUE, then the data is centered following the \code{\link{prcomp}} function. If it is set to TRUE, scale input is ignored.
  }
\item{scale}{
If it is set to TRUE, the data is scaled through the spectral mean value. Set to FALSE by default.
  }
\item{signVariables}{
If this input field is a numeric vector, only the spectra/peaks having an ID number present in such vector are used as input data. If it is set to NULL, all the variables are taken into account to build the input data.
  }
\item{RemoveOnePeakSpectra}{
If it is set to TRUE, all the one-peak spectra are deleted from the dataSet and the resulting spectralData object will only contain spectra with more than one peak.
  }
\item{printCSVfile}{
If it is set to TRUE, an output matrix showing the spectral/peak intensity is build, where each column is a sample and each row is a variable (spectra or peak depending on the method used).
  }

}
\value{
  An \code{\link{MAIT-class}} object.
}
\examples{
data(MAIT_sample)
peakAggregation(MAIT)
}

\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}
\keyword{function}
\keyword{internal}