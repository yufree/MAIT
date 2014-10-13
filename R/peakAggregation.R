########################################################################################################################################################## 
##########################################################################################################################################################
##                                                      |
##   Metabolite Automatic Identification Toolkit (MAIT) |
##                                                      |
##                                                      |
##   written by Francesc Fernández Albert               |
##   contact mail: francesc.fernandez.albert@upc.edu    |
##   date: 7/29/2013                                    |
##                                                      |
##   SISBIO Group. ESAII Department                     |
##   Technical University of Catalonia                  |
##   Nutrition Department                               |
##   University of Barcelona                            |
##   ___________________________________________________|
##
##    
##   peakAggregation function applies one of five possible peak aggregation techniques (none, spectral mean, single peak, PCA decomposition and NMF reduction.
##   The needed variables to run the function are:
##
##   MAIT.object:          MAIT object where the data is taken from. The output object of the function is going to be an update of the same MAIT object.
##   method:               Method used for dimensionality reduction of the data.
##   clases:               Classes present in the xsaFA variable.
##   samples:              If it is a numeric vector, only the samples with the ID number present in the vector are chosen (set to NULL by default).
##   PCAscale:             If TRUE scale of the function prcomp is used instead of the scale done by getDataSet function (FALSE by default).
##   PCAcenter:            If TRUE center of the function prcomp is used instead of the scale done by getDataSet function (FALSE by default).
##   scale:                Data is scaled dividing each variable by its mean value.
##   signVariables:        If it is a numeric vector, only the variables with the ID number present in the vector are chosen (set to NULL by default).
##   RemoveOnePeakSpectra: If it is TRUE, the spectra with just one peak are removed from the dataSet (set to FALSE by default).
##   printCSVfile:         If it is set to TRUE, a CSV file is printed under the name of (working directory)/Tables/dataSet.csv.
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
##########################################################################################################################################################
##########################################################################################################################################################
peakAggregation<-function(MAIT.object=NULL,
                          method="None",
                          clases=NULL,
                          samples=NULL,
                          PCAscale=FALSE,
                          PCAcenter=FALSE,
                          scale=FALSE,
                          signVariables=NULL,
                          RemoveOnePeakSpectra=FALSE,
                          printCSVfile=TRUE
                          ){
if (is.null(MAIT.object)){
        stop("No input MAIT object was given")
    }

if (is.null(method)){
        stop("No input peak aggregation method was given")
    }

if(method=="NMF"|method=="PCA"|method=="Mean"|method=="Single"){
        stop("Your MAIT version do not support using peak aggregation measures. Make sure that pagR package is installed and that you are using the correct MAIT version.")
      }

        parameters <- list(method,
                           PCAscale,
                           PCAcenter,
                           scale,
                           signVariables,
                           RemoveOnePeakSpectra)

        names(parameters) <- c("peakAggregation method",
                           "peakAggregation PCAscale",
                           "peakAggregation PCAcenter",
                           "peakAggregation scale",
                           "peakAggregation signVariables",
                           "peakAggregation RemoveOnePeakSpectra")
       MAIT.object@RawData@parameters@peakAggregation <- parameters
       writeParameterTable(parameters(MAIT.object),folder=resultsPath(MAIT.object))
  classes <- classes(MAIT.object)
  classNum <- classNum(MAIT.object)
  resultsPath <- resultsPath(MAIT.object)

  MAIT.object@FeatureInfo@peakAgMethod <- method
##########################################################################################################################################################
#
#   Peaks intensity is gathered from the xsaFA object and ordered by its spectra id number. idGroup stores the spectra id information
#
##########################################################################################################################################################

aux <- getScoresTable(MAIT.object=MAIT.object,getSpectra=TRUE,getExtendedTable=FALSE)
   data <- aux$scores
    if(scale==TRUE && method!="Single"){
      
       data<-data/rowMeans(data)
       data[is.nan(as.matrix(data))]<-0

}
idGroup <- aux$spectraID 

##########################################################################################################################################################
#
#   If RemoveOnePeakSpectra is set to TRUE function removeOnePeakSpectra is invoked on the data, removing all the one-peak spectra.
#
##########################################################################################################################################################
        if(RemoveOnePeakSpectra==TRUE){
           dataWithNoOnePeakSpectra <- removeOnePeakSpectra(data=data,
                                                     idGroup=idGroup)
           data <- dataWithNoOnePeakSpectra$spectra
           idGroup <- dataWithNoOnePeakSpectra$idGroup
         }
##########################################################################################################################################################
#
#   If samples is a numeric vector, just the samples whose id number are in the vector are selected.
#
##########################################################################################################################################################
       if(!is.null(samples)){
          data <- data[,samples]
       }
##########################################################################################################################################################
#
#   Nas coming from the xcms peak dectection step are set to 0.
#
##########################################################################################################################################################

    data[is.na(as.matrix(data))]<-0
    colnames(data) <- NULL
    rownames(data) <- NULL
##########################################################################################################################################################
#
#   If signVariables is set to NULL, all the variables are used to perform the peak aggregation and build the data set.
#
##########################################################################################################################################################
 if(is.null(signVariables)){

   MAIT.object@FeatureData@scores <- data
   MAIT.object@FeatureData@featureID <- idGroup

     MAIT.object@FeatureData@scores <- as.matrix(data)
     MAIT.object@FeatureData@featureID <- idGroup
     MAIT.object@FeatureData@models <- list(rep(1,dim(data)[1]))

     
##########################################################################################################################################################
#
#   If signVariables is a numeric vector, only the variables whose id number are in the vector are selected in the data set.
#    
##########################################################################################################################################################
  }else{
    index <- vector(length=1)
    noneIdGroup <- vector(length=1)
    for(i in c(1:length(signVariables))){
        index<-c(index,which(signVariables[i]==idGroup))
        noneIdGroup<-c(noneIdGroup,rep(signVariables[i],length(which(signVariables[i]==idGroup))))
    }
 index <- sort(index)
 index <- index[-1]
 noneIdGroup <- noneIdGroup[-1]


     MAIT.object@FeatureData@scores <- as.matrix(data[index,])
     MAIT.object@FeatureData@featureID <- noneIdGroup
     MAIT.object@FeatureData@models <- list(rep(1,dim(data)[1]))
     MAIT.object@FeatureInfo@peakAgMethod <- "None"
  }
##########################################################################################################################################################
#
#   In the following, the so-called dataSet.csv table is printed if the input parameter printCSVfile is set to TRUE
#    
##########################################################################################################################################################
  if(printCSVfile==TRUE){
       if(!file.exists(paste(resultsPath,"Tables",sep="/"))){

         dir.create(paste(resultsPath,"Tables",sep="/"))
       }else{
         cat(" " ,fill=TRUE)
#         warning(paste("Folder",paste(resultsPath,"Tables",sep="/"),"already exists. Possible file overwritting.",sep=" "),fill=TRUE)
       }

       tabl <- scores(MAIT.object)
       if(length(rawData(MAIT.object))!=0){
  if(is.null(samples)){
       colnames(tabl) <- sampnames(rawData(MAIT.object)[[1]]@xcmsSet)

  }else{
    colnames(tabl) <- sampnames(rawData(MAIT.object)[[1]]@xcmsSet)[samples]
  }
 
}else{
  colnames(tabl) <- colnames(scores(MAIT.object))
}
             rownames(tabl) <- paste("S",1:dim(tabl)[1],sep="")
       if (scale==TRUE){
         norm <- "scaled"
       }else{
         norm <- "notScaled"
       }
       write.table(file=paste(resultsPath,paste("Tables/dataSet_",norm,".csv",sep=""),sep="/"),x=tabl,row.names=TRUE,col.names=NA,sep=",")
    }

    return(MAIT.object)

  }
