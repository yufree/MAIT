MAITbuilder<-function(data=NULL,
		      spectraID=NULL,
                      masses=NULL,
                      rt=NULL,
                      classes=NULL,
                      significantFeatures=FALSE,
                      spectraEstimation=FALSE,
                      rtRange=0.2,
                      corThresh=0.7){


  
  
  if(is.null(data) & length(masses)==0){

stop("No peak data or masses set as an input")

}
  if(is.null(data)){
    warning("No input peak data was given")
  }

  
MAIT.object<-new("MAIT")
  if(is.null(data)!=TRUE){
if(is.null(colnames(data))){
colnames(data) <- paste("Sample",1:dim(data)[2])
}
}
  if(is.null(classes)==TRUE|is.null(classNum)==TRUE){
    MAIT.object@PhenoData@classes<-NA
    MAIT.object@PhenoData@classNum<-NA
  }else{
    classNum <- summary(as.factor(classes[order(classes)]))
    data <- data[,order(classes)]
    classes <- unique(classes[order(classes)])
    MAIT.object@PhenoData@classes <- classes    
    MAIT.object@PhenoData@classNum <- classNum
  }

     if(is.null(masses)!=TRUE){
MAIT.object@FeatureData@masses <- masses
}else{
    if(is.null(data)!=TRUE){
masses <- rep(NA,nrow(data))
}else{
  masses <- NA 
}
}
   if(is.null(rt)!=TRUE){
MAIT.object@FeatureData@rt <- rt
}else{
    if(is.null(data)!=TRUE){
rt <- rep(NA,nrow(data))
}else{
  rt <- rep(NA,length(masses))
}
}

   if(is.null(spectraID)==TRUE){
     if(is.null(data)!=TRUE){
         if(spectraEstimation==FALSE){
spectraID <- 1:dim(data)[1]
}else{
    spectraID <- computeSpectra(peaks=data, rt=rt,rtRange=rtRange, corThresh=corThresh)
  spectraID <- as.numeric(factor(spectraID))}

}else{
  if(spectraEstimation==FALSE){
  spectraID <- 1:length(masses)
}else{
  spectraID <- computeSpectra(peaks=data, rt=rt,rtRange=rtRange, corThresh=corThresh)
  spectraID <- as.numeric(factor(spectraID))
}
}
}
  
  if(is.null(data)!=TRUE){
MAIT.object@FeatureData@scores <- data
}else{
data <- rep(NA,length(masses))
}






  aux <- unique(data.frame(data,masses,rt,spectraID))


data <- subset(aux,select=-c(masses,spectraID,rt))
  rt <- aux$rt
  masses <- aux$masses
  spectraID <- aux$spectraID
MAIT.object@FeatureData@scores <- as.matrix(data)
MAIT.object@FeatureData@rt <- rt
MAIT.object@FeatureData@masses <- masses

  

tab <- as.data.frame(cbind(matrix(masses,ncol=1),matrix(rt,ncol=1),data,spectraID))
  colnames(tab) <- c("mz","rt",colnames(data))
  colnames(tab)[dim(tab)[2]] <- "pcgroup"
  tab <- tab[order(tab$pcgroup),]
MAIT.object@FeatureData@extendedTable <- tab

  if(significantFeatures==TRUE){
    MAIT.object@FeatureData@featureSigID <- tab$pcgroup

  }

  
return(MAIT.object)
}
