########################################################################################################################################################## 
##########################################################################################################################################################
##                                                      |
##   Metabolite Automatic Identification Toolkit (MAIT) |
##                                                      |
##                                                      |
##   written by Francesc Fernández Albert               |
##   contact mail: francesc.fernandez.albert@upc.edu    |
##   date: 10/15/2012                                   |
##                                                      |
##   SISBIO Group. ESAII Department                     |
##   Universitat Politècnica de Catalunya               |
##   ___________________________________________________|
##
##    
##   Validation function performs a repeated random sub-sampling cross-validation step using three different classifiers: KNN, PLSDA and SVM.
##
##   Iterations:           Number of classifications to perform
##   trainSamples:         Number of samples per class used to construct the train dataset
##   class:                Class names involved in the xsAnnotate object
##   classNum:             Number of elements belonging to each class
##   MAIT.object:          Object containing the spectral data after having applied the function spectralAnova or spectralTStudent over it.
##   PCAscale:             Scale of the function prcomp is used instead of the scale done by getDataSet function (FALSE by default)
##   PCAcenter:            Center of the function prcomp is used instead of the scale done by getDataSet function (FALSE by default)
##   signData:             Output of the function signPeaks signData
##   method:               Chosen method to reduce the dimensionality of the data ("NMF" by default)
##   RemoveOnePeakSpectra: Spectra with only one peak are removed from the dataSet (FALSE by default)
##
##
## All code copyright (c) 2013 UPC/UB
## All accompanying written materials copyright (c) 2013 UPC/UB
##
##
##   This program is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
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


Validation <- function(Iterations=NULL,
                           MAIT.object=NULL,
                           trainSamples=NULL,
                           PCAscale=FALSE,
                           PCAcenter=TRUE,
                           RemoveOnePeakSpectra=FALSE,
                           tuneSVM=FALSE,
                           scale=TRUE){
 
    if (is.null(Iterations)) {
        stop("No Iterations number set!")
    }
  
    if (is.null(MAIT.object)) {
        stop("No input MAIT object file was given")
    }

     if (is.null(trainSamples)) {
       stop("No input trainSamples parameter was given.")
    }

     if (length(classes(MAIT.object))==1) {
       if(is.na(classes(MAIT.object))==TRUE){
       stop("No class information in the MAIT object.")
    }
     }






    
        parameters <- list(Iterations,
                           trainSamples,
                           PCAscale,
                           PCAcenter,
                           RemoveOnePeakSpectra,
                           tuneSVM,
                           scale)

            names(parameters) <- c("Validation Iterations",
                           "Validation trainSamples",
                           "Validation PCAscale",
                           "Validation PCAcenter",
                           "Validation RemoveOnePeakSpectra",
                           "Validation tuneSVM",
                           "Validation scale")


       MAIT.object@RawData@parameters@classification <- parameters
        writeParameterTable(parameters(MAIT.object),folder=resultsPath(MAIT.object))


    
method <- method(MAIT.object)
class<- classes(MAIT.object)
classNum <- classNum(MAIT.object)
resultsPath <- resultsPath(MAIT.object)

##########################################################################################################################################################
#
#    Variable definition. SRKnn, SRPLSDA and SRSVM will hold the classification ratios.
#
##########################################################################################################################################################


        SRKnn <- matrix(nrow=Iterations)
	SRPLSDA <- matrix(nrow=Iterations)
	SRSVM <- matrix(nrow=Iterations)
        SR_class <- matrix(nrow=Iterations,ncol=3*length(class))
  
    
	groupTrain <- vector(length=1)
	group <- vector(length=1)
    
	for(i in c(1:length(class))){
		groupTrain <- c(groupTrain,rep(class[i],trainSamples))
		group <- c(group,rep(class[i],classNum[i]))
	}
    
	groupTrain <- as.factor(groupTrain[-1])
	group <- as.factor(group[-1])


    
	classAc <- c(1,classNum[1])
	for (i in c(2:length(classNum))){
		classAc <- c(classAc,classAc[i]+classNum[i])
	}
  
        
##########################################################################################################################################################
#        
#     ClassWeights takes into account the population number differences between classes. This variable is used to modify the classification ratio.
#
##########################################################################################################################################################

        
	ClassWeights <- vector(length=length(class))
	for (i in c(1:length(classNum))){
		ClassWeights[i] <- 1-(classNum[i])/sum(classNum)
          	
	}

  

##########################################################################################################################################################
#
#     Loop over the desired Iterations number
#
##########################################################################################################################################################

        
         k<-1
    
         if(!file.exists(paste(resultsPath,"Validation",sep="/"))){
								   
            dir.create(paste(resultsPath,"Validation",sep="/"))
            dir.create(paste(resultsPath,"Validation","Confusion_Tables",sep="/"))
                                								   
          }else{
								   
	    warning(paste("Folder",paste(resultsPath,"Validation",sep="/"),"already exists. Possible file overwritting.",sep=" "))
								   
	  }

    
          while(k<=Iterations){

          
          
##########################################################################################################################################################
#
#  Definition of the training and validation indexs through a random choice
#
##########################################################################################################################################################

          
		trainIndex <- matrix(nrow=1)
                trainIndex <- c(trainIndex,sample(classAc[2]:(classAc[1]),trainSamples,replace=FALSE))
		for (i in c(3:length(classAc))){
			trainIndex <- c(trainIndex,sample(classAc[i]:(classAc[i-1]+1),trainSamples,replace=FALSE))
		}
		trainIndex <- trainIndex[-1]

                valIndex<-1:classAc[length(classAc)]
                for (i in c(1:length(trainIndex))){
                valIndex<-valIndex[-which(valIndex==trainIndex[i])]
                
              }
 
              
              



                
##########################################################################################################################################################
#
#  Dimensionality reduction is performed while the method chosen is not "None". Variable trainDataSet is a spectralData object whereas train contains the scores of the
#  training set.
#
##########################################################################################################################################################
                


                   

                       peaks <- getScoresTable(MAIT.object)$scores
                
                       peaks[is.na(as.matrix(peaks))]<-0
                       train <- peaks[featureSigID(MAIT.object),trainIndex]
                       train<-train/apply(train,1,mean)
                       train[is.nan(as.matrix(train))]<-0
                       val <- peaks[featureSigID(MAIT.object),valIndex]
                       val<-val/apply(val,1,mean)
                       val[is.nan(as.matrix(val))]<-0

                

##########################################################################################################################################################
#
#    Three diferent classifications are KNN, PLSDA and SVM are performed in the workflow
#                
##########################################################################################################################################################
              
#                
#
                
##########################################################################################################################################################
#                
#    PLSDA Classifier
#
##########################################################################################################################################################
                
		PLScomps <- try(selectPLScomp(data=t(train),class=groupTrain,max.comp=trainSamples*length(class)))
                if(!is.numeric(PLScomps)){PLScomps <- 1}
                if(trainSamples>2){
                model <- try(plsda(t(train),groupTrain,ncomp=PLScomps,scale=TRUE),TRUE)
              }else{
                model <- try(plsda(t(train),groupTrain,ncomp=PLScomps,scale=FALSE),TRUE)}

                if(!is.character(model[1])){
                                  
		   pred <- predict(object=model,newdata=t(as.data.frame(val)))
		   tt <- table(pred=pred,true=group[-trainIndex])
                   

                   writeExcelTable(tt,paste(resultsPath,"Validation","Confusion_Tables",paste("PLSDAConfTable_",k,sep=""),sep="/"))

		   succRatioPLSDA <- successRatio(classes=class,tt=tt,ClassWeights=ClassWeights)
                   SRPLSDA[k] <- succRatioPLSDA$SR
                 }
                   
##########################################################################################################################################################
#                
#    KNN Classifier
#
##########################################################################################################################################################
                 
                   
                    
                     
		Kcomps <- selectK(t(train),class=groupTrain,max.k=(trainSamples*length(class))-1)
		knnRaw <- knn(train=t(train),test=t(val),cl=groupTrain,k=Kcomps)
		tt <- table(pred=knnRaw,true=group[-trainIndex])

                writeExcelTable(tt,paste(resultsPath,"Validation","Confusion_Tables",paste("KNNConfTable_Iteration_",k,sep=""),sep="/"))
                   
		succRatioKNN <- successRatio(class,tt,ClassWeights=ClassWeights)
                SRKnn[k] <- succRatioKNN$SR

              
                   
              	Kcomps <- selectK(t(train),class=groupTrain,max.k=(trainSamples*length(class))-1)
		knnRaw <- knn(train=t(train),test=t(val),cl=groupTrain,k=Kcomps)
		tt <- table(pred=knnRaw,true=group[-trainIndex])

                writeExcelTable(tt,paste(resultsPath,"Validation","Confusion_Tables",paste("KNNConfTable_",k,sep=""),sep="/"))
                   
		succRatioKNN <- successRatio(class,tt,ClassWeights=ClassWeights)
                SRKnn[k] <- succRatioKNN$SR



                
##########################################################################################################################################################
#                
#    SVM Classifier
#
##########################################################################################################################################################

                
                if(tuneSVM==TRUE){
                  
                   tun<-tune.svm(t(train),groupTrain, gamma = 10^(-6:-3), cost = 10^(1:3))
                   if (trainSamples>1){
		   modelSVM <- (svm(t(train),groupTrain,type="C",kernel="radial",scale=TRUE,gamma=tun$best.parameters$gamma,cost=tun$best.parameters$cost))
                 }else{
                   modelSVM <- (svm(t(train),groupTrain,type="C",kernel="radial",scale=FALSE,gamma=tun$best.parameters$gamma,cost=tun$best.parameters$cost))
                 }

                   
                }else{
                   if (trainSamples>1){
		   modelSVM <- svm(t(train),groupTrain,type="C",kernel="radial",scale=TRUE,gamma=0.001,cost=10)
                 }else{
                   modelSVM <- svm(t(train),groupTrain,type="C",kernel="radial",scale=FALSE,gamma=0.001,cost=10)
                                 }
                 }
                   
		   pred <- predict(modelSVM,t(val))
		   tt <- table(pred=pred,true=group[-trainIndex])
                   writeExcelTable(tt,paste(resultsPath,"Validation","Confusion_Tables",paste("SVMConfTable_",k,sep=""),sep="/"))
		   succRatioSVM <- successRatio(classes=class,tt=tt,ClassWeights=ClassWeights)
                   
                   SRSVM[k] <- succRatioSVM$SR

              
                   
                    cat(paste("Iteration",as.character(k),"done",sep=" "),fill=TRUE);
                   
                       if (k<=Iterations){

                           
                          for(i in c(1:length(class))){

                             SR_class[k,(i*3)-2] <- succRatioKNN$SR_class[i] 
                             SR_class[k,(i*3)-1] <- succRatioPLSDA$SR_class[i]
                             SR_class[k,(i*3)] <- succRatioSVM$SR_class[i]
                          
                        
                       
                       }
                 
                   k<-k+1
                       }
              }
                   
  
                
            
##########################################################################################################################################################
#                
#    Output is generated and saved in a folder called (working directory)/Validation. Mean and standard error are computed from the classification iterations. 
#
##########################################################################################################################################################

                
    	       out <- matrix(ncol=3,nrow=2)
               names_class <- vector(length=3*length(class))
               out_class <- matrix(ncol=3,nrow=2)
               colour <- vector(length=3*length(class))


for (i in c(1:length(class))){

  
               out_class[1,1] <- mean(SR_class[,(i*3)-2])
    	       out_class[2,1] <- sd(SR_class[,(i*3)-2])/sqrt(ncol(train))
	       out_class[1,2] <- mean(SR_class[,(i*3)-1])
   	       out_class[2,2] <- sd(SR_class[,(i*3)-1])/sqrt(ncol(train))
	       out_class[1,3] <-mean(SR_class[,(i*3)])
	       out_class[2,3] <- sd(SR_class[,(i*3)])/sqrt(ncol(train))
    	       colnames(out_class) <- c("KNN","PLSDA","SVM")
               rownames(out_class) <- c("mean","standard error")
  	       writeExcelTable(out_class,paste(resultsPath,"Validation",paste("ClassificationTable_Class_",class[i],sep=""),sep="/"))

               names_class[(i*3)-2] <- paste("KNN_Class_",class[i],sep="")
               names_class[(i*3)-1] <- paste("PLSDA_Class_",class[i],sep="")
               names_class[(i*3)] <- paste("SVM_Class_",class[i],sep="")
               
}

  
               palette <- palette(rainbow(length(class)))
    
               for (i in c(1:length(palette))){
                  colour[((i*3)-2):(i*3)] <- rep(palette[i],3)
               }

               colnames(SR_class) <- rep(c("KNN","PLSDA","SVM"),length(class))

               
               png(paste(resultsPath,"Validation","Boxplot_Clases_Classification.png",sep="/"))
               boxplot(SR_class,horizontal=TRUE,col=colour,pars=list(box=FALSE,axis=FALSE),xlab="Weighted classification ratio",las=1)
               legend("topright",legend=class, pt.bg=palette,bty='n',pch=22)
               title("Classification ratio per class depending on the classifier")
               dev.off()

               colnames(SR_class) <- names_class


    
           #    out[1,1] <- mean(SRKnn)
    	    #   out[2,1] <- sd(SRKnn)/sqrt(ncol(train))
	    #   out[1,2] <- mean(SRPLSDA)
   	    #   out[2,2] <- sd(SRPLSDA)/sqrt(ncol(train))
	    #   out[1,3] <- mean(SRSVM)
	    #   out[2,3] <- sd(SRSVM)/sqrt(ncol(train))

               out[1,1] <- mean(SRKnn)
    	       out[2,1] <- apply(SRKnn,2,sd)/sqrt(ncol(train))
	       out[1,2] <- mean(SRPLSDA)
   	       out[2,2] <- apply(SRPLSDA,2,sd)/sqrt(ncol(train))
	       out[1,3] <- mean(SRSVM)
	       out[2,3] <- apply(SRSVM,2,sd)/sqrt(ncol(train))

    	       colnames(out) <- c("KNN","PLSDA","SVM")

               rownames(out) <- c("mean","standard error")

                
##########################################################################################################################################################
#
#   A boxplot and some tables showing the results (overall and for each class) are created in the following
#
##########################################################################################################################################################

                
  	       writeExcelTable(out,paste(resultsPath,"Validation","ClassificationTable",sep="/"))
    
	       outBox <- list(SRKnn,SRPLSDA,SRSVM)
	       names(outBox) <- c("KNN","PLSDA","SVM")

	       png(paste(resultsPath,"Validation","Boxplot_Overall_Classification.png",sep="/"))
	       boxplot(outBox,ylab="Weighted Classification Accuracy",horizontal=TRUE)
               title("Overall classification ratio")

	       dev.off()

##########################################################################################################################################################
#
#   The results of the search are saved in the provided MAIT object which is the output of the function.
#
##########################################################################################################################################################

                
               MAIT.object@Validation@classifRatioClasses <- SR_class
               MAIT.object@Validation@ovClassifRatioTable <- out
               MAIT.object@Validation@ovClassifRatioData <- outBox
	       return(MAIT.object)
             }
