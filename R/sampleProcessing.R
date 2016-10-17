########################################################################################################################################################## 
##########################################################################################################################################################
##                                                      |
##   Metabolite Automatic Identification Toolkit (MAIT) |
##                                                      |
##                                                      |
##   written by Francesc Fernández Albert               |
##   contact mail: francesc.fernandez.albert@upc.edu    |
##   date: 7/29/2013                                   |
##                                                      |
##   SISBIO Group. ESAII Department                     |
##   Technical University of Catalonia                  |
##   Nutrition Department                               |
##   University of Barcelona                            |
##                                                      |
##   SISBIO Group. ESAII Department                     |
##   Universitat Politècnica de Catalunya               |
##   ___________________________________________________|
##
##
##   sampleProcessing function takes a set of netCDF files containing LC/MS sample data and performs a peak detection, retention time correction and peak           ##   grouping steps using the xcms package. A MAIT object is created and all the informated is saved in it. Each of these steps are done using the package xcms.    ##   The input variables of the function are:
##
##   dataDir:       The netCDF sample files of each class present in the data should be stored in a folder called /(working directory)/Data/(ClassName) replacing   ##                  (ClassName) for the name of the folder where the files are stored.
##   snThres:       Signal to noise ratio. Setting a high value of this parameter will lead to a higher number of features although they will be more noisy         ##                  (set to 2 by default).
##   Sigma:         Standard deviation (width) of matched filtration model peak.
##   mzSlices:      Minimum difference in m/z for peaks with overlapping retention times.
##   retcorrMethod: Method used to correct the retention times values of the variables. By default is set to "loess".
##   groupMethod:   Method used to build the group peaks of variables. By default is set to "density".
##   bwGroup:       Bandwidth (standard deviation or half width at half maximum) of gaussian smoothing kernel to apply to the peak density chromatogram.
##   mzWidGroup:    Width of overlapping m/z slices to use for creating peak density chromatograms and grouping peaks across samples.
##   filterMethod:  Filtering method applied in the peak detection step. (Set to "matchedFilter" by default).
##   rtStep:        Step size to use for profile generation.
##   nSlaves:       Number of slaves for parallel calculus.
##   project:       Project folder name under which the results will be saved. This folder will be created in the working directory.
##
##
## All code copyright (c) 2013 UPC/UB
## All accompanying written materials copyright (c) 2013 UPC/UB
##
##
##   This program is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 2 of the License, or
##    (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with this program.  If not, see <http://www.gnu.org/licenses/>.
##
##    
########################################################################################################################################################## 
##########################################################################################################################################################


sampleProcessing <- function(dataDir=NULL,
                             snThres=2,
                             Sigma=5/2.3548,
                             mzSlices=0.3,
                             retcorrMethod="loess",
                             groupMethod="density",
                             bwGroup=3,
                             mzWidGroup=0.25,
                             filterMethod="matchedFilter",
                             rtStep=0.03,
                             nSlaves=0,
                             minfrac=0.5,
                             minsamp=1,
                             peakwidth=c(5,20),
			     project=NULL,
                             ppm=10,
                             family=c("gaussian", "symmetric"),
                             span = .2,
                             fwhm=30)
{

   

  if (is.null(dataDir)) {
        stop("No input directory was given")
    }

    if (is.null(project)) {
        stop("No project name was included")
    }



       parameters <- list(dataDir,
                             snThres,
                             Sigma,
                             mzSlices,
                             retcorrMethod,
                             groupMethod,
                             bwGroup,
                             mzWidGroup,
                             filterMethod,
                             rtStep,
                             nSlaves,
			     project,
                             ppm,
                             minfrac,
                             fwhm,
                             family,
                             span,
                             peakwidth)
  
       names(parameters) <- c("dataDir",
                             "snThres",
                             "Sigma",
                             "mzSlices",
                             "retcorrMethod",
                             "groupMethod",
                             "bwGroup",
                             "mzWidGroup",
                             "filterMethod",
                             "rtStep",
                             "nSlaves",
			     "project",
                             "ppm",
                             "minfrac",
                             "fwhm",
                             "family",
                              "span",
                              "centWave peakwidth")
  
       MAIT.object <- new("MAIT")

  
       MAIT.object@RawData@parameters@sampleProcessing <- parameters
       writeParameterTable(parameters(MAIT.object),folder=resultsPath(MAIT.object))
  
       class <-  list.files(dataDir)
       classNum <- vector(length=length(class)) 
       fileList <- list.files(path=paste(dataDir,list.files(path=dataDir),sep="/"),full.names=TRUE)

       for (i in 1:length(class)){
          classNum[i] <- length(list.files(paste(dataDir,list.files(dataDir)[i],sep="/")))
       }

  
       classes <- rep(class,classNum)

    if (length(list.files(dataDir))==1){
    warning("Warning: Input data only has one class!")
  }
	
	if (is.null(project)){
		warning("Warning: Project name is empty!")
	}

##########################################################################################################################################################
#
#  The results folder is created if it does not exist
#
##########################################################################################################################################################


  
if (!is.null(project)) {
		
		resultsPath<-paste("Results",project,sep="_")
		dir.create(resultsPath)

	}else{
		
		resultsPath<-"Results"
		dir.create(resultsPath)
		
	}	

  

        
##########################################################################################################################################################
#
#  Peak detection step.
#
##########################################################################################################################################################

        if (filterMethod=="matchedFilter"){
	peaks<- xcmsSet(files=fileList,snthresh=snThres,method=filterMethod,sigma=Sigma,max=3,step=rtStep,mzdiff=mzSlices,sclass=classes,nSlaves=nSlaves,fwhm=fwhm)
      }
        if (filterMethod=="centWave"){
        peaks<- xcmsSet(files=fileList,snthresh=snThres,method=filterMethod,ppm=ppm,mzdiff=mzSlices,sclass=classes,nSlaves=nSlaves,peakwidth=peakwidth)
      }
	cat("Peak detection done",fill=TRUE)

        
##########################################################################################################################################################
#
#  Peak grouping and retention time correction steps.
#
##########################################################################################################################################################

        
	groups<- group(peaks,method=groupMethod,bw=bwGroup,mzwid=mzWidGroup,max=50,minfrac=minfrac,minsamp=minsamp)
  
  if(retcorrMethod!="none"){
    
	retcorr_groups<- retcor(groups,method=retcorrMethod,plottype="deviation",family = family,span=span)
	cat("Retention time correction done",fill=TRUE)
	groups<- group(retcorr_groups,method=groupMethod,bw=bwGroup,mzwid=mzWidGroup,max=50)
	cat("Peak grouping after samples done",fill=TRUE)
        
  }else{
    
        cat("Skipping retention time correction...",fill=TRUE)

  }

        
##########################################################################################################################################################
#
#  Fill missing peaks step
#
##########################################################################################################################################################

        
	fPeaks<- fillPeaks(groups)
	cat("Missing Peak integration done",fill=TRUE)

##########################################################################################################################################################
#
#  Writting the results in a new MAIT object
#
##########################################################################################################################################################
        fPeaks <- list(fPeaks)
        names(fPeaks) <- "xcmsSet"

  
        MAIT.object@RawData@data <- fPeaks        
        MAIT.object@PhenoData@classes <- class
        MAIT.object@PhenoData@classNum <- classNum
        MAIT.object@PhenoData@resultsPath <- resultsPath

        return(MAIT.object)
}

