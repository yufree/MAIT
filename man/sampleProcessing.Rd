\name{sampleProcessing}
\alias{sampleProcessing}
\title{Peak detector of netCDF samples using xcms package}
\description{
  sampleProcessing takes a set of netCDF files containing LC/MS sample data and performs a peak detection, retention time correction and peak grouping steps using the package xcms. A \link{MAIT-class} object is created and all the informated is saved in it.
}
\usage{
  sampleProcessing(dataDir = NULL, 
                   snThres = 2, 
                   Sigma = 5/2.3548,
		   mzSlices = 0.3, 
		   retcorrMethod = "loess", 
		   groupMethod = "density", 
		   bwGroup = 3, 
		   mzWidGroup = 0.25,
		   filterMethod = "matchedFilter",
  		   rtStep = 0.03, 
		   nSlaves = 0, 
		   minfrac = 0.5,
                   minsamp = 1, 
		   peakwidth = c(5, 20), 
		   project = NULL, 
		   ppm = 10,
                   family = c("gaussian", "symmetric"),
                   span = 0.2,
                   fwhm = 30)
}
\arguments{
  \item{dataDir}{
Folder where the netCDF files are stored. The samples files must be classified in subdirectories according to their classes.
}
  \item{snThres}{
Signal to noise ratio. Setting a high value of this parameter will lead to a higher number of features although they will be more noisy.
}
  \item{Sigma}{
    Standard deviation (width) of matched filtration model peak.
  }
  \item{mzSlices}{
     Minimum difference in m/z for peaks with overlapping retention times.
  }
  \item{retcorrMethod}{
     Method used to correct the retention times values of the variables.
  }
\item{groupMethod}{
     Method used to build the group peaks of variables.
  }
\item{bwGroup}{
    Bandwidth (standard deviation or half width at half maximum) of gaussian smoothing kernel to apply to the peak density chromatogram.
  }
\item{mzWidGroup}{
    Width of overlapping m/z slices to use for creating peak density chromatograms and grouping peaks across samples.
  }
\item{filterMethod}{
    Filtering method applied in the peak detection step.
  }
\item{rtStep}{
   Step size to use for profile generation.
  }
\item{nSlaves}{
   Number of slaves for parallel calculus.
  }
\item{project}{
   Project folder name under which the results will be saved. This folder will be created in the working directory.
  }
  \item{minfrac}{
    minimum fraction of samples necessary in at least one of the sample groups for it to be a valid group. See group.density in package \pkg{xcms} for details.
  }
  \item{minsamp}{
    minimum number of samples necessary in at least one of the sample groups for it to be a valid group. See group.density in package \pkg{xcms} for details.
  }
\item{ppm}{
	maxmial tolerated m/z deviation in consecutive scans, in ppm (parts per million). See findPeaks.centWave in package \pkg{xcms} for details.
	}
\item{peakwidth}{
	Chromatographic peak width, given as range (min,max) in seconds.
	}
  \item{fwhm}{
    See fwhm argument in xcmsSet function.
  }

  \item{span}{
    See span argument in xcmsSet function.
  }
  \item{family}{
    See family argument in xcmsSet function.
  }	
}
\value{
  A \link{MAIT-class} object containing the data of the netCDF files. The xcmsSet-class object can be retrieved using the function \link{rawData}.
}

\examples{
#Provided that the data files are saved accordingly 
#in subfolders under a folder named "data" (see vignette):
#MAIT<-sampleProcessing(dataDir = "data", project = "Results", snThres=2,rtStep=0.02)
}

\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}
\keyword{function}
