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

  
sigPeaksTable<-function(
		MAIT.object=NULL,
                printCSVfile=FALSE,
                extendedTable=TRUE,
                printAnnotation=TRUE){
 

 
##########################################################################################################################################################
##       
##   sigPeaksTable is an output table build of some features of statistically significant variables. This table is build from a MAIT object.
##
########################################################################################################################################################## 


   if (is.null(MAIT.object)) {
        stop("No MAIT object was given")
    }


if(length(featureSigID(MAIT.object))==0){
          stop("No significant features found in the MAIT object. Make sure that functions peakAnnotation and spectralSigFeatures were launched")
        }

   
#        peakList <- getPeaklist(MAIT.object)
 #  	peakList <- peakList[order(as.numeric(peakList$pcgroup)),]

        dat <- getScoresTable(MAIT.object,getSpectra=TRUE,getExtendedTable=TRUE)

   if(printAnnotation==TRUE){extendedTable<-TRUE}
   
   if(extendedTable==TRUE){
        peakList <- dat$extendedTable
      }else{
        peakList <- dat$scores
      }
        spectraID <- dat$spectraID
   
        data <- scores(MAIT.object)
        index <- featureSigID(MAIT.object)
        classes <- classes(MAIT.object)
        classNum <- classNum(MAIT.object)
        resultsPath <- resultsPath(MAIT.object)
        Fisher <- LSDResults(MAIT.object)
        TTs <- pvalues(MAIT.object)

   

	sigPeaksTable <- matrix(nrow=1,ncol=ncol(peakList))
	colnames(sigPeaksTable) <- colnames(peakList)

        
	if (length(index)!=0){
          

           if (method(MAIT.object)!="None"){
             
              sigPeaksTable <- peakList[spectraID%in%index,]
             
           }else{
          
             sigPeaksTable <- peakList[spectraID%in%unique(spectraID[index]),]
             
           }

           
           
	   p <- matrix(nrow=1,ncol=dim(sigPeaksTable)[1])
	   fisher <- matrix(nrow=1,ncol=dim(sigPeaksTable)[1])

                      if(class(classes(MAIT.object))!="logical"|class(classNum(MAIT.object))!="logical"){

           
           if(length(classes(MAIT.object))>2){
           
	   for(i in c(1:dim(sigPeaksTable)[1])){
          
              fisher[i] <- Fisher[as.numeric(sigPeaksTable[i,dim(sigPeaksTable)[2]])]
           
	   }

           	   fisher <- as.vector(fisher)
	   fisherNames <- paste(classes[1],classes[2],sep="_")
           
	   for (i in c(3:length(classNum))){
          
	      fisherNames <- paste(fisherNames,classes[i],sep="_")
        
	   }
           
	   NamesFisher <- paste("Fisher",as.character(fisherNames),sep=" ")
	   names(fisher) <- NamesFisher
           

           
         }else{
           
           for(i in c(1:dim(sigPeaksTable)[1])){
          
              fisher[i] <- NA
           
	   }
         }

         }
           
           p<-pvalues(MAIT.object)
           p <- p[spectraID%in%unique(spectraID[index])]

	   p <- matrix(p,ncol=1)
           if(MAIT.object@FeatureData@pvaluesCorrection==""){MAIT.object@FeatureData@pvaluesCorrection<-"none"}
	   P.adjust <- p.adjust(p,MAIT.object@FeatureData@pvaluesCorrection)

           pAux <- matrix(ncol=1,nrow=dim(sigPeaksTable)[1])
           pAux.adjust <- matrix(ncol=1,nrow=dim(sigPeaksTable)[1])
 

           if (method(MAIT.object)!="None"){
             
           for (i in c(1:dim(sigPeaksTable)[1])){
             if(sum(sigPeaksTable$pcgroup==index)==length(index)){
           pAux[i,] <- unique(p[which(index%in%index[i])])
           pAux.adjust[i,] <- unique(P.adjust[which(index%in%index[i])])

         }else{
                      pAux[i,] <- p[which(index%in%as.numeric(sigPeaksTable$pcgroup)[i])]
           pAux.adjust[i,] <- P.adjust[which(index%in%as.numeric(sigPeaksTable$pcgroup)[i])]
                    }
           }

         }else{

           pAux <- p
           pAux.adjust <- P.adjust

         }
         
	   sigPeaksTable <- cbind(sigPeaksTable,pAux.adjust,pAux)

#   if(length(classes(MAIT.object))>2){

     sigPeaksTable <- cbind(sigPeaksTable,matrix(fisher,ncol=1))
              if(length(classes(MAIT.object))>2){
     colnames(sigPeaksTable)[dim(sigPeaksTable)[2]] <- NamesFisher
   }else{
     colnames(sigPeaksTable)[dim(sigPeaksTable)[2]] <- "Fisher.Test"
   }
     	   colnames(sigPeaksTable)[dim(sigPeaksTable)[2]-2] <- "P.adjust"
       	   colnames(sigPeaksTable)[dim(sigPeaksTable)[2]-1] <- "p"
 #  }else{
     
#	   colnames(sigPeaksTable)[dim(sigPeaksTable)[2]-1] <- "P.adjust"
 #      	   colnames(sigPeaksTable)[dim(sigPeaksTable)[2]] <- "p"
  #  }
 

           

if(length(rawData(MAIT.object))==0&is.null(sigPeaksTable$adduct)==TRUE){
  adduct <- character(length=dim(peakList)[1])
  peakList <- data.frame(peakList,adduct)
    adduct <- character(length=dim(sigPeaksTable)[1])

  sigPeaksTable <- data.frame(sigPeaksTable,adduct)
sigPeaksTable$adduct<-as.character(sigPeaksTable$adduct)
  peakList$adduct<-as.character(peakList$adduct)
}

  
  
if(sum(is.na(MAIT.object@FeatureInfo@biotransformations))!=1&dim(MAIT.object@FeatureInfo@biotransformations)[1]!=0){
  
aux<-as.data.frame(MAIT.object@FeatureInfo@biotransformations)
aux[,1]<-as.numeric(unlist(MAIT.object@FeatureInfo@biotransformations[,1]))
aux[,2]<-as.numeric(unlist(MAIT.object@FeatureInfo@biotransformations[,2]))
  
biotrans <- aux

for(i in c(1:dim(biotrans)[1])){
  if(sigPeaksTable[as.numeric(biotrans[i,7]),]$adduct==""){

    sigPeaksTable[as.numeric(biotrans[i,7]),]$adduct <- as.character(biotrans[i,3])
    
  }else{
    
sigPeaksTable[as.numeric(biotrans[i,7]),]$adduct <- paste(sigPeaksTable[as.numeric(biotrans[i,7]),]$adduct,biotrans[i,3],sep=";")

}
}
}

   
   
sigPeaksTable$adduct <- unlist(sigPeaksTable$adduct)

           means <- vector("list",length=length(index))
           medians <- vector("list",length=length(index))

           if(class(classes(MAIT.object))!="logical"|class(classNum(MAIT.object))!="logical"){
           
           Fgroups <- as.factor(rep(classes(MAIT.object),classNum(MAIT.object)))
             
           for(i in c(1:dim(sigPeaksTable)[1])){

             if(length(rawData(MAIT.object))==0){
               ind <- c(3:(dim(sigPeaksTable)[2]-5))
             }else{
               ind <- c((8+length(classes(MAIT.object))):(dim(sigPeaksTable)[2]-6))
             }
               
       #     means[[i]] <- aggregate(dat$scores[spectraID%in%unique(spectraID[index]),][i,]~Fgroups,FUN=mean)
       #     medians[[i]] <- aggregate(dat$scores[spectraID%in%unique(spectraID[index]),][i,]~Fgroups,FUN=median)
             means[[i]] <- aggregate(as.numeric(sigPeaksTable[i,ind])~Fgroups,FUN=mean)
             medians[[i]] <- aggregate(as.numeric(sigPeaksTable[i,ind])~Fgroups,FUN=median)
          }

           allMeans <- merge(means[[1]],means[[2]],by="Fgroups")
           if(length(rawData(MAIT.object))!=0){
           colnames(allMeans)[2:dim(allMeans)[2]] <- rownames(dat$scores[spectraID%in%unique(spectraID[index]),])[c(1,2)]
         }else{
#           colnames(allMeans)[2:dim(allMeans)[2]] <- rep(NA,length(2:dim(allMeans)[2]))
           colnames(allMeans)[2:dim(allMeans)[2]] <- 1:(dim(allMeans)[2]-1)
         }
           if(length(index)>2){
           for(i in c(3:length(means))){
             allMeans <- merge(allMeans,means[[i]],by="Fgroups")
             if(length(rawData(MAIT.object))!=0){
                colnames(allMeans)[i+1] <-  rownames(dat$scores[spectraID%in%unique(spectraID[index]),])[i]
              }else{
                colnames(allMeans)[i+1]  <- as.numeric(colnames(allMeans)[i])+1
           }
           }
         }
            # colnames(allMeans)[2:dim(allMeans)[2]] <- rownames(dat$scores)
             rownames(allMeans)<-paste("Mean Class",allMeans[,1])

           temp<-as.data.frame(t(allMeans[,-1]))
           for(i in c(1:length(classes(MAIT.object)))){
           temp[,i]<-as.numeric(as.character(temp[,i]))
         }

          if(length(rawData(MAIT.object))!=0){rownames(temp)<-NULL}

           allMedians <- merge(medians[[1]],medians[[2]],by="Fgroups")
          if(length(rawData(MAIT.object))!=0){
           colnames(allMedians)[2:dim(allMedians)[2]] <- rownames(dat$scores[spectraID%in%unique(spectraID[index]),])[c(1,2)]
                    }else{
#           colnames(allMeans)[2:dim(allMeans)[2]] <- rep(NA,length(2:dim(allMeans)[2]))
           colnames(allMedians)[2:dim(allMedians)[2]] <- 1:(dim(allMedians)[2]-1)
         }
           if(length(index)>3){
           for(i in c(3:length(medians))){
             allMedians <- merge(allMedians,medians[[i]],by="Fgroups")
             if(length(rawData(MAIT.object))!=0){
             colnames(allMedians)[i+1] <-  rownames(dat$scores[spectraID%in%unique(spectraID[index]),])[i]
                         }else{
                colnames(allMedians)[i+1]  <- as.numeric(colnames(allMedians)[i])+1

           }
           }
           }
            # colnames(allMedians)[2:dim(allMedians)[2]] <- rownames(dat$scores)
             rownames(allMedians)<-paste("Median Class",allMedians[,1])

           tempMed<-as.data.frame(t(allMedians[,-1]))
           for(i in c(1:length(classes(MAIT.object)))){
           tempMed[,i]<-as.numeric(as.character(tempMed[,i]))
         }
           sigPeaksTable <- cbind(sigPeaksTable,temp,tempMed)

          if(length(rawData(MAIT.object))!=0){rownames(tempMed)<-NULL}


         }else{

           temp <- matrix(rep(NA,dim(sigPeaksTable)[1]),ncol=1)
           tempMed <- matrix(rep(NA,dim(sigPeaksTable)[1]),ncol=1)

                      sigPeaksTable <- cbind(sigPeaksTable,temp,tempMed)
colnames(sigPeaksTable)[c(dim(sigPeaksTable)[2]-1,dim(sigPeaksTable)[2])] <- c("Class Mean", "Class Median")
           
         }
           

           
           if(printCSVfile==TRUE){

				if(!file.exists(paste(resultsPath,"Tables",sep="/"))){
					
					dir.create(paste(resultsPath,"Tables",sep="/"))

				}else{

					cat(" " ,fill=TRUE)
					warning(paste("Folder",paste(resultsPath,"Tables",sep="/"),"already exists. Possible file overwritting.",sep=" "))

				}


	           write.csv(x=sigPeaksTable,file=paste(resultsPath,"Tables/significantFeatures.csv",sep="/"),row.names=FALSE)

            }

           
       	}else{
             
           stop("There are no significant features for this pvalue")
           
         }

   
 return(sigPeaksTable)
 }
 
