computeSpectra <- function(peaks=NULL, rt=NULL, rtRange=NULL, corThresh=NULL){

if(is.null(peaks)|is.null(rt)){

stop("Not enough data was provided!")

}

corPeaks<-cor(t(peaks))
diffrt <- diff(rt)

spectraID <- 1:length(rt)

i <- 1
spectraID[1]=1
index <- 1:length(rt)
for(i in c(1:length(spectraID))){
  if(!length(which(spectraID[i]==spectraID))>1){
    if(length(which(abs(rt[i]-rt[c(which(corPeaks[i,]>=corThresh)[-i])])<rtRange))>0){
      ind <- which(corPeaks[i,]>=corThresh)[-i][which(abs(rt[i]-rt[c(which(corPeaks[i,]>=corThresh)[-i])])<rtRange)]
      spectraID[ind] <- spectraID[i]
      index <- index[-ind]
    }
  }
}

return(spectraID)
}
