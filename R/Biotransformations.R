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
##   Biotransformations function looks for specifical relationships between the peaks of the same spectra. These relationships are provided by the user through a   ##   csv table.
##
##   MAIT.object:          MAIT object where the data is taken from. The output object of the function is going to be an update of the same MAIT object.
##   peakPrecision:        In order to find the peak relationships, this parameter computes the tolerance of the peak mass differences. As a consequence,           ##                         the candidate interval for each pair of peaks is going to be [peak1-peak2-peakPrecision,peak1-peak2+peakPrecision].
##                         All the entries of the bioTable csv file falling within these margins will be considered as a suitable biotransformation.
##   bioTable:             This parameter must be a csv table having two columns. Each entry (row) must be a possible biotransformation whereas the first column    ##                         should be the name of the biotransformation and the second column, the mass difference of such biotransformation.
##
##
##   All code copyright (c) 2013 UPC/UB 
##   All accompanying written materials copyright (c) 2013 UPC/UB
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

 
Biotransformations <- function(MAIT.object=NULL,
                               peakPrecision=0.005,
                               bioTable=NULL,
                               adductTable=NULL,
                               adductAnnotation=FALSE){
     
    if (is.null(MAIT.object)) {
        stop("No input MAIT object file was given")
    }

      if(length(featureSigID(MAIT.object))==0){
          stop("No significant features found in the MAIT object. Make sure that functions peakAnnotation and peakAggregation were launched")
        }

        if(length(MAIT.object@FeatureData@masses)==0 & length(rawData(MAIT.object))==0){
      stop("No peak masses found in the MAIT object")
    }
    
          if (is.null(bioTable)) {
              peakBioEnv<-new.env()
              data(MAITtables,envir=peakBioEnv)
              biotransformationsTable<-get("biotransformationsTable",envir=peakBioEnv)
              
              Biotransformations <- biotransformationsTable
#        Biotransformations <- biotransformationsTable
              warning("No input biotransformations table was given. Selecting default MAIT table for biotransformations...")
        }else{
         biotransformationsTable <- as.data.frame(read.csv(paste(bioTable,".csv",sep=""),sep=",",header=TRUE,row.names=NULL))
       	cat("Selecting user-defined Biotransformations Table",fill=TRUE)

      }
       
H.mass <- 1.00794
if(adductAnnotation){

    if (is.null(adductTable)){

      	    warning("No input adduct/fragment table was given. Selecting default MAIT table for positive polarity...")
            cat("Set adductTable equal to negAdducts to use the default MAIT table for negative polarity",fill=TRUE)
            peakAnnEnv<-new.env()
            data(MAITtables,envir=peakAnnEnv)
            adducts<-get(x="posAdducts",envir=peakAnnEnv)
            adducts$massdiff<-adducts$massdiff-adducts$charge*H.mass
     }else{

    if (adductTable=="negAdducts"){
    cat("adductTable has been set to negAdducts. The default MAIT adducts table for negative polarization is selected...",fill=TRUE)

                peakAnnEnv<-new.env()
                data(MAITtables,envir=peakAnnEnv)
      adducts<-get(x="negAdducts",envir=peakAnnEnv)
      adducts$massdiff<-adducts$massdiff-adducts$charge*H.mass
    
    }else{

         adducts <- read.csv2(paste(adductTable,".csv",sep=""),dec=".",header=TRUE,sep=",")
         
       }
  }
    
    
    temp <- adducts[which(adducts$nmol==1),c("name","massdiff")]
    colnames(temp) <- colnames(biotransformationsTable)
    biotransformationsTable <- rbind(biotransformationsTable,temp)

adducts$name <- as.character(adducts$name)
}
    
   Biotransformations <- biotransformationsTable

       parameters <- list(peakPrecision,
                          bioTable,
                          adductTable,
                          adductAnnotation)

       names(parameters) <- c("peakPrecision",
                              "bioTable",
                              "adductTable",
                              "Biotransformations adductAnnotation")


       MAIT.object@RawData@parameters@biotransformations <- parameters
       writeParameterTable(parameters(MAIT.object),folder=resultsPath(MAIT.object))
    
##########################################################################################################################################################
#
#   Getting the required data from the MAIT object.
#
##########################################################################################################################################################


    if(length(rawData(MAIT.object))!=0){
    
        if(!names(rawData(MAIT.object))=="xsaFA"){
          
          stop("No peak annotation has been performed yet")
          
        }
      }
            sigPeaksTable <- sigPeaksTable(MAIT.object)

	BiotransformationsTable <- matrix(ncol=7)
	aux <- as.numeric(biotransformationsTable[,2])
	index <- matrix(nrow=1)
    	colnames(BiotransformationsTable) <- c("Mass1","Mass2","Adduct","Biotransf Name","rt","Initial Peak","Final Peak")

            
        cat('\n % Annotation in progress: ')
        lp <- -1

##########################################################################################################################################################
#
#   signIndex corresponds to the significant spectra IDs
#
##########################################################################################################################################################

    
signIndex <- sort(unique(as.numeric(sigPeaksTable$pcgroup)))


  
biotRange <- cbind(aux-peakPrecision,aux+peakPrecision)
rownames(biotRange)<-biotransformationsTable$NAME
colnames(biotRange) <- c("Margin-","Margin+")

#if(adductPolarity=="positive"){
#addRange <- cbind(adducts$massdiff-1-peakPrecision,adducts$massdiff-1+peakPrecision,adducts$nmol)
#}else{
#  addRange <- cbind(adducts$massdiff+1-peakPrecision,adducts$massdiff+1+peakPrecision,adducts$nmol)
#}

    
allSpectra <- sapply(signIndex,FUN=retrieveSpectrum,sigPeaksTable=sigPeaksTable)


    
if(adductAnnotation){
    
#for(i in c(1:length(as.numeric(levels(factor(adducts$nmol)))))){
    addRange <- cbind(adducts$massdiff-peakPrecision,adducts$massdiff+peakPrecision,adducts$nmol)
    
    colnames(addRange) <- c("Margin-","Margin+","nmol")
    rownames(addRange) <- adducts$name

    auxTab<-adducts[which(adducts$nmol==1),c("name","massdiff")]
#    auxTab<-adducts[,c("name","massdiff")]
#colnames(auxTab) <- c("NAME","MASSDIFF","NMOL")

for(k in 1:length(signIndex)){

  spectrum <- unique(unlist(allSpectra[[k]]))
if(length(spectrum)>1){
names(spectrum)<-match(round(spectrum,4),round(sigPeaksTable$mz,4))
#out <- addDiff(x=spectrum,nmol=as.numeric(levels(factor(adducts$nmol)))[i])
#diffList <- lapply(c(1:(length(spectrum)-1)),FUN=diff,x=spectrum)
BiotransformationsTable <- annotAdduct(spectrum=spectrum,addRange=addRange,BiotransformationsTable=BiotransformationsTable,auxTab=auxTab,sigPeaksTable=sigPeaksTable)
diffList <- lapply(c(1:(length(spectrum)-1)),FUN=diff,x=spectrum)
#BiotransformationsTable <- spectrumAnnotation(testValues=diffList,BiotransformationsTable=BiotransformationsTable,spectrum=spectrum,biotRange=addRange[which(addRange[,3]==1),],sigPeaksTable=sigPeaksTable,biotransformationsTable=auxTab)
#addRange[which(adducts$nmol==i),]
#adducts[which(adducts$nmol==i),]
}
}
}
    





    
addDiff <- function(x,nmol){
out <- vector("list",length=length(x)-1)
  aux <- x
for(i in c(1:(length(x)-1))){
  aux[i] <- aux[i]*nmol
  out[[i]] <- lapply(c(1:(length(aux)-1)),FUN=diff,x=aux)
  aux[i] <- x[i]
   }
return(out)
}






















##########################################################################################################################################################
#
#   Starting the loop for all the significant features
#
##########################################################################################################################################################
 
    
for(k in 1:length(signIndex)){
spectrum <- unique(unlist(allSpectra[[k]]))
if(length(spectrum)>1){
names(spectrum)<-match(round(spectrum,4),round(sigPeaksTable$mz,4))
diffList <- lapply(c(1:(length(spectrum)-1)),FUN=diff,x=spectrum)

BiotransformationsTable <- spectrumAnnotation(testValues=diffList,BiotransformationsTable=BiotransformationsTable,spectrum=spectrum,biotRange=biotRange,sigPeaksTable=sigPeaksTable,biotransformationsTable=biotransformationsTable)
}
spectrum<-rev(unlist(spectrum))
if(length(spectrum)>1){
names(spectrum)<-match(round(spectrum,4),round(sigPeaksTable$mz,4))
diffList <- lapply(c(1:(length(spectrum)-1)),FUN=diff,x=spectrum)
BiotransformationsTable <- spectrumAnnotation(testValues=diffList,BiotransformationsTable=BiotransformationsTable,spectrum=spectrum,biotRange=biotRange,sigPeaksTable=sigPeaksTable,biotransformationsTable=biotransformationsTable)
}
perc <-round(which(signIndex[k]==signIndex)/length(signIndex)*100)
if ((perc %% 10 == 0) && (perc != lp)) { cat(perc,' '); lp <- perc }
}
BiotransformationsTable <- (BiotransformationsTable[-1,])
if(class(BiotransformationsTable)!="list"){
        MAIT.object@FeatureInfo@biotransformations <- BiotransformationsTable
      }else{
        aux<-names(BiotransformationsTable)
        BiotransformationsTable <- matrix(BiotransformationsTable,nrow=1,ncol=7)
        colnames(BiotransformationsTable) <- aux
        MAIT.object@FeatureInfo@biotransformations <- BiotransformationsTable
      }

  
    	return(MAIT.object)
  }

            
##########################################################################################################################################################
#
#   Function spectrumAnnotation performs the whole Biotransformation annotation stage for a single spectrum
#
##########################################################################################################################################################



spectrumAnnotation <- function(testValues,BiotransformationsTable,spectrum,biotRange,sigPeaksTable=sigPeaksTable,biotransformationsTable,diffSpectrum=NULL){
out <- vector("list",length=length(testValues))
for(i in c(1:length(testValues))){

  if(is.null(dim(sapply(X=testValues[[i]],FUN=inBetween,biotRange=biotRange)))){
    out[[i]] <- as.list(sapply(X=testValues[[i]],FUN=inBetween,biotRange=biotRange))
  }else{
    auxNames <- vector()
    for(j in c(1:dim(sapply(X=testValues[[i]],FUN=inBetween,biotRange=biotRange))[2])){
      auxNames <- c(auxNames,rep(colnames(as.list(sapply(X=testValues[[i]],FUN=inBetween,biotRange=biotRange)))[j],dim(sapply(X=testValues[[i]],FUN=inBetween,biotRange=biotRange))[1]))
     #      auxNames <- c(auxNames,rep(colnames(aux)[j],dim(sapply(X=testValues[[i]],FUN=inBetween,biotRange=biotRange))[1]))
    }
  }

if(!is.null(diffSpectrum)){testValues <- diffSpectrum}
  
if(length(unlist(out[[i]]))!=0){
  counter <- 0
  indices <- integer()
for(o in c(1:length(out[[i]]))){
  if(!identical(unlist(out[[i]][o]), character(0))){
    counter <- counter+1;
    indices <- c(indices,o)
  }
}
  for(j in c(1:counter)){
BiotransformationsTable <- rbind(BiotransformationsTable,annotateBiotransf(biotransf=(out[[i]])[indices[j]],diffIndex=i,spectrum=spectrum,sigPeaksTable=sigPeaksTable,biotransformationsTable=biotransformationsTable))
}
}
}
return(BiotransformationsTable)
}

##########################################################################################################################################################
#
#   Function inBetween looks up if a mass difference might correspond to a certain Biotransformation
#
##########################################################################################################################################################


inBetween<-function(testValue,biotRange){
  out <- rownames(biotRange)[which(testValue<=biotRange[,2]&testValue>=biotRange[,1])]
  names(out) <- rep(names(testValue),length(out))
  return(out)
}

##########################################################################################################################################################
#
#   Function annotateBiotransf annotates a single biotransformation for a certain spectrum
#
##########################################################################################################################################################

  
 annotateBiotransf <- function(biotransf,diffIndex,spectrum,sigPeaksTable,biotransformationsTable){
  if(length(unlist(biotransf))>1){
  biotransfMass <- biotransformationsTable$MASS[match(table=biotransformationsTable$NAME,x=unlist(biotransf)[1])]
  biotransfName <- character(length=1)
  for(i in c(1:length(unlist(biotransf)))){
    biotransfName <- paste(biotransfName,unlist(biotransf)[i],",")
  }
}else{
  biotransfMass <- biotransformationsTable$MASS[match(table=biotransformationsTable$NAME,x=biotransf)]
  }
  auxIndex2 <- as.numeric(names(spectrum)[which(names(spectrum)%in%names(biotransf))])
  auxIndex1 <- as.numeric(names(spectrum)[which(names(spectrum)%in%names(biotransf))-diffIndex])
  massPeak2 <- sigPeaksTable$mz[auxIndex2]
  massPeak1 <- sigPeaksTable$mz[auxIndex1]
  rtPeak2 <- sigPeaksTable$rt[match(table=sigPeaksTable$mz,x=sigPeaksTable$mz[auxIndex2])]
  signPeakIndex2 <- match(table=sigPeaksTable$mz,x=sigPeaksTable$mz[auxIndex2])
  rtPeak1 <- sigPeaksTable$rt[match(table=sigPeaksTable$mz,x=sigPeaksTable$mz[auxIndex1])]
  signPeakIndex1 <- match(table=sigPeaksTable$mz,x=sigPeaksTable$mz[auxIndex1])
  if(biotransfMass>0){
   if(length(unlist(biotransf))>1){
     annotation <- paste("[M",round(as.numeric(massPeak1),3),"_rt",round(as.numeric(rtPeak1),2),"] ",biotransfName,sep="")
     output <- matrix(c(massPeak1,massPeak2,annotation,biotransfName,rtPeak1,signPeakIndex1,signPeakIndex2),nrow=1,ncol=7)
   }else{
     annotation <- paste("[M",round(as.numeric(massPeak1),3),"_rt",round(as.numeric(rtPeak1),2),"] ",as.character(biotransf),sep="")
     output <- matrix(c(massPeak1,massPeak2,annotation,biotransf,rtPeak1,signPeakIndex1,signPeakIndex2),nrow=1,ncol=7)
   }

  }else{
    if(length(unlist(biotransf))>1){
     annotation <- paste("[M",round(as.numeric(massPeak1),3),"_rt",round(as.numeric(rtPeak1),2),"] ",biotransfName,sep="")
     output <- matrix(c(massPeak1,massPeak2,annotation,biotransfName,rtPeak2,signPeakIndex1,signPeakIndex2),nrow=1,ncol=7)
    }else{
     annotation <- paste("[M",round(as.numeric(massPeak1),3),"_rt",round(as.numeric(rtPeak1),2),"] ",as.character(biotransf),sep="")
     output <- matrix(c(massPeak1,massPeak2,annotation,biotransf,rtPeak2,signPeakIndex1,signPeakIndex2),nrow=1,ncol=7)

        }
  }
  return(output)
}

##########################################################################################################################################################
#
#   Function retrieveSpectrum retrieves the peak list for a certain spectrum ID
#
##########################################################################################################################################################


 retrieveSpectrum <- function(spectrumNumber,sigPeaksTable){return(sigPeaksTable[which(as.numeric(sigPeaksTable$pcgroup)==spectrumNumber),1])}

annotAdduct <- function(spectrum,addRange,BiotransformationsTable,auxTab,sigPeaksTable){

if(max(unique(addRange[,3]))>1){
for(i in c(2:max(unique(addRange[,3])))){
for(j in c(1:length(spectrum))){
  addRangeMod <- addRange[which(addRange[,3]==i),]
  aux <- abs(spectrum[j]*i-spectrum[-j])
names(aux) <- names(spectrum)[-j]
  if(!is.null((sapply(X=aux,FUN=inBetween,biotRange=addRangeMod)))){
    an <- sapply(X=aux,FUN=inBetween,biotRange=addRangeMod)
    for(k in c(1:length(an))){
    if(length(unlist(an[k]))>=1){
      mass1 <- sigPeaksTable$mz[as.numeric(names(an))[k]]
      mass2 <- sigPeaksTable$mz[as.numeric(names(spectrum))[j]]
      rtPeak1 <- sigPeaksTable$rt[as.numeric(names(an))[k]]
      rtPeak2 <- sigPeaksTable$rt[as.numeric(names(spectrum))[j]]
      signPeakIndex2 <- match(table=sigPeaksTable$mz,x=mass1)
      signPeakIndex1 <- match(table=sigPeaksTable$mz,x=mass2)
 
#     annotation <- paste("[M",round(as.numeric(massPeak1),3),"_rt",round(as.numeric(rtPeak1),2),"] ",biotransfName,sep="")
      annotation <- paste(an[k]," M",round(mass2,4),"_rt",rtPeak1,sep="")
      biotransfName <- an[k]
     output <- matrix(c(mass1,mass2,annotation,biotransfName,rtPeak1,signPeakIndex1,signPeakIndex2),nrow=1,ncol=7)

      BiotransformationsTable <- rbind(BiotransformationsTable,output)
    }
  }
  }
}
}
}
return(BiotransformationsTable)
}
