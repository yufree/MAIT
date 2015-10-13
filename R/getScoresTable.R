getScoresTable <- function(MAIT.object=NULL,
                           getSpectra=TRUE,
                           getExtendedTable=FALSE){

  
   if (is.null(MAIT.object)) {
        stop("No MAIT object was given")
    }

if(length(rawData(MAIT.object))==0){

peakList <- MAIT.object@FeatureData@extendedTable
scores <- scores(MAIT.object)
spectra <- MAIT.object@FeatureData@extendedTable$pcgroup
     out <- list(scores,as.numeric(spectra),peakList)
     names(out) <- c("scores","spectraID","extendedTable")

}else{
   
   classes<- classes(MAIT.object)

   peakList <- getPeaklist(MAIT.object)
   peakList<-peakList[order(as.numeric(peakList$pcgroup)),]

   if(length(which(classes[length(classes)]==colnames(peakList))+1)!=0){
   scores <- as.matrix(peakList[,(which(classes[length(classes)]==colnames(peakList))+1):(length(colnames(peakList))-3)])
}else{
   scores <- as.matrix(peakList[,(which(paste("X",classes[length(classes)],sep="")==colnames(peakList))+1):(length(colnames(peakList))-3)])
 }
   
   spectra <- peakList$pcgroup

     if(getExtendedTable==TRUE){

     if(getSpectra==TRUE){

           peakList[,4:6] <- round(peakList[,4:6]/60,2)
       
     out <- list(scores,as.numeric(spectra),peakList)
     names(out) <- c("scores","spectraID","extendedTable")

    }else{

     out <- list(scores,peakList)
     names(out) <- c("scores","cameraTable")
     
   }
   }else{
     
          if(getSpectra==TRUE){
       
     out <- list(scores,as.numeric(spectra))
     names(out) <- c("scores","spectraID")

    }else{

     out <- list(scores)
     names(out) <- c("scores")
     
   }        
   }
 }
   
   return(out)

 }

