\name{Validation}
\alias{Validation}
\title{Cross validated classification over the output of the function \code{signPeaksAnova} or function \code{TStudent2Clases}}
\description{
Function Validation performs a cross-validated classification using three different classifiers: KNN, PLSDA and SVM. The output comes in a table with the classification ratio and its standard error. The classification ratio is weighted to take into account the different sample number of each class.
}
\usage{

Validation(Iterations=NULL,
               MAIT.object=NULL,
               trainSamples=NULL,
               PCAscale=FALSE,
               PCAcenter=TRUE,
               RemoveOnePeakSpectra=FALSE,
               tuneSVM=FALSE,
	           scale=TRUE)
}

\arguments{
  \item{Iterations}{
Number of iterations to be performed in the classifications. For each iteration a new training group is randomly chosen.
}
  \item{MAIT.object}{
A \link{MAIT-class} object where significant features have already been found.
}
  \item{trainSamples}{
Number of samples per class to construct the train dataset.
}
  \item{PCAscale}{
If method="PCA" and PCAscale is set to TRUE, then the data is scaled following the \code{prcomp} function. If it is set to TRUE, scale input is ignored.
  }
\item{PCAcenter}{
If method="PCA" and PCAscale is set to TRUE, then the data is centered following the \code{prcomp} function. If it is set to TRUE, scale input is ignored.
  }
\item{RemoveOnePeakSpectra}{
If it is set to TRUE, all the one-peak spectra are deleted from the dataSet and the resulting \code{spectralData} object will only contain spectra with more than one peak.
  }
\item{tuneSVM}{
If it is set to TRUE, a tune of parameters is performed before the SVM calculus.
}
\item{scale}{
If it is set to TRUE, the data is scaled through the spectral mean value. Set to TRUE by default.
}


}
\value{
The numerical results of the classification per class and per classifier are saved in a \link{MAIT-class} object. Additionally, a table is also included in the output both in the list (field table) and printed as a csv file in the folder (working directory)/Validation. A boxplot is also printed as a png in the same folder showing the differences between classifiers. The confusion matrices of each iteration and classifier are also stored as csv files.
}
\seealso{
  \code{\link{peakAggregation}}  \code{\link{spectralAnova}} \code{\link{spectralTStudent}} \code{\link{spectralSigFeatures}} }


\examples{
data(MAIT_sample)
MAIT<-spectralSigFeatures(MAIT,p.adj="fdr",parametric=TRUE)
MAIT <- Validation(Iterations = 20, trainSamples= 15, MAIT.object = MAIT)
}


\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}
