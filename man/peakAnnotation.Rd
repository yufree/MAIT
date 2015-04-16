\name{peakAnnotation}
\alias{peakAnnotation}
\title{Spectra constructor and peak annotator}
\description{
  peakAnnotation function performs spectra building and peak annotation using the CAMERA package on a \link{MAIT-class} object, after applying the \link{sampleProcessing} function. The resultant \link{xsAnnotate} object is stored in a \link{MAIT-class} object.
}
\usage{
peakAnnotation(MAIT.object = NULL, 
               corrWithSamp = 0.7, 
               perfwhm = 0.6, 
               sigma = 6, 
               adductTable = NULL, 
               printSpectraTable = TRUE, 
               corrBetSamp = 0.75,
               pval = 0.05, 
               calcIso = TRUE,
               calcCiS = TRUE, 
               calcCaS = TRUE,
               graphMethod = "hcs", 
               annotateAdducts = TRUE)
}

\arguments{
  \item{MAIT.object}{
A \link{MAIT-class} object where function \link{sampleProcessing} has already been applied. The output of the function is going to be an update of the same \link{MAIT-class} object.
}
  \item{corrWithSamp}{
Correlation threshold value within samples
}
  \item{perfwhm}{
This parameter is used to group two peaks depending on their retention time. Two peaks are considered to be coeluted if their retention time falls in a range defined as  Rt_med +/- FWHM * perfwhm. Where Rt_med is the retention time median and FWHM is the Full Width at Half Maximum. Defined this way, perfwhm is the percentage of the width of the FWHM (Full Width at Half Maximum)}
  \item{sigma}{
Defining the coelution range as defined in the perfwhm variable, the FWHM is obtained by the expression FWHM=SD*sigma, where SD is calculated considering the peak as normally distributed.
 }
  \item{adductTable}{
User-defined input table to annotate the peaks. If it is set to NULL, the default MAIT table for adducts in positive polarization is selected. If its value is "negAdducts", the default MAIT table for fragments in negative polarization is chosen. By default it is set to NULL. 
  }
\item{printSpectraTable}{
If it is set to TRUE, a three-column table is build as a csv file, where the first column shows the peak mass, the second column its retention time and the third one shows its spectral ID number. This file is saved under the project directory, in the subfolder named Tables.
  }
  \item{corrBetSamp}{
Correlation threshold value between samples
}
  \item{pval}{
See \link{groupCorr} function in the \pkg{CAMERA} package
}
  \item{calcIso}{
See \link{groupCorr} function in the \pkg{CAMERA} package
}
  \item{calcCiS}{
See \link{groupCorr} function in the \pkg{CAMERA} package
}
  \item{calcCaS}{
See \link{groupCorr} function in the \pkg{CAMERA} package
}
  \item{graphMethod}{
See \link{groupCorr} function in the \pkg{CAMERA} package
}
  \item{annotateAdducts}{
If it is set to TRUE, the function will perform an adduct annotation stage.
}


}
\value{
  A \link{MAIT-class} object containing the \link{xsAnnotate-class} in the rawData slot.
}
\seealso{
  \code{\link{xsAnnotate},\link{xsAnnotate-class}}
}

\examples{
#Provided that the data files are saved accordingly 
#in subfolders under a folder named "data" (see vignette):
#MAIT<-sampleProcessing(dataDir = "data", project = "Results", snThres=2,rtStep=0.02)
#MAIT<-peakAnnotation(MAIT.object = MAIT,corrWithSamp = 0.7, corrBetSamp = 0.7,perfwhm = 0.6)
}

\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}
\keyword{file}
