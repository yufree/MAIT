\name{MAITbuilder}
\alias{MAITbuilder}
\title{MAIT constructor function when using external peak data}
\description{
Function MAITbuilder creates a MAIT-class object for a given external data. This process allows the user to analyse external peak data through all the MAIT processing steps.
}
\usage{

MAITbuilder(data = NULL, 
            spectraID = NULL,
            masses = NULL, 
            rt = NULL, 
            classes = NULL, 
            significantFeatures = FALSE, 
            spectraEstimation = FALSE, 
            rtRange = 0.2,
            corThresh = 0.7)
}

\arguments{
  \item{data}{
Matrix containing the peak intensity values for each sample. Each row should correspond to a peak and each column to a sample.
}
  \item{spectraID}{
Numeric corresponding to the peak spectral grouping IDs.Two peaks having the same spectraID means that they correspond to the same spectrum.
}
  \item{masses}{
Numeric that contains the masses of the peaks. It should be as long as the number of rows in the argument data.
}
  \item{rt}{
Numeric that contains the retention time of the peaks. It should be as long as the number of rows in the argument data.
  }
\item{classes}{
Character with the class labels for each sample. It should be as long as the number of columns in the argument data.
  }
\item{significantFeatures}{
If it is set to TRUE, all the features set as an input are considered to be significant. Funcions \link{Biotransformations}, \link{identifyMetabolites}, \link{Validation}, \link{plotPCA}, \link{plotPLS}, \link{plotHeatmap}, \link{plotBoxplot} are computed on the significant features only. If it is only wanted to perform an annotation process on the external peak data, this flag should be set to TRUE.
  }
\item{spectraEstimation}{
If it is set to TRUE, an estimation of the peak grouping into spectra is performed. This computation is based on a retention time window (set by the argument rtRange) and a correlation threshold (defined by the parameter corThresh).
}
\item{rtRange}{
Retention time parameter used to build a window to perform an estimation of the peak grouping into spectra.
}
\item{corThresh}{
Peak correlation value used to define a threshhol to perform an estimation of the peak grouping into spectra.
}

}

\value{
All the imput values are stored in a new MAIT object.
}

\examples{
data(MAIT_sample)
peaks<-scores(MAIT)
aux<-getScoresTable(MAIT)
masses<-aux$extendedTable$mz
rt <- aux$extendedTable$rt
classFactor <- rep(classes(MAIT),classNum(MAIT))
importMAIT <- MAITbuilder (data=peaks,masses=masses,rt=rt,
significantFeatures=TRUE, spectraEstimation=TRUE, rtRange=0.2,
corThresh=0.7,classes=classFactor)

importMAIT
}


\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}
