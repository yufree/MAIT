########################################################################################################################################################## 
##########################################################################################################################################################
##                                                      |
##   Metabolite Automatic Identification Toolkit (MAIT) |
##                                                      |
##                                                      |
##   written by Francesc Fernández Albert               |
##   contact mail: francesc.fernandez.albert@upc.edu    |
##   date: 1/20/2014                                    |
##                                                      |
##   SISBIO Group. ESAII Department                     |
##   Technical University of Catalonia                  |
##   Nutrition Department                               |
##   University of Barcelona                            |
##   ___________________________________________________|
##  
##      
##   Function identifyMetabolites takes the data from a MAIT object
##   defined to build an output table made of significant features only. Fisher LSD analysis is also performed in order to know which are the class differences in the
##   statistically significant features.
##   
##   MAIT.object:    Object containing the spectral processed data (coming from function getDataSet). ANOVA tests will be performed on each of its features.
##   peakTolerance:  In order to find the peak relationships, this parameter computes the tolerance of the peak mass differences. As a consequence, the candidate
##                   interval for each pair of peaks is going to be [peak1-peak2-peakPrecision,peak1-peak2+peakPrecision]. All the entries of the bioTable csv file falling within these margins will be
##                   considered as a suitable biotransformation.
##   bonferroni:     If it is set to TRUE, Bonferroni multiple testing correction is performed. Pvalue threshold is taken after the correction
##   printCSVfile:   If it is set to TRUE, a table containing the results of the statistical tests is printed under the name
##                   (working directory)/Tables/significativeFeatures.csv through the sigPeaksTable function
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
  
 
identifyMetabolites <- function(MAIT.object=NULL,
                                peakTolerance=0.005,                                
                                database=NULL,
                                polarity="positive",
                                printCSVfile=TRUE){
  


    if (is.null(MAIT.object)) {
        stop("No input MAIT object file was given")
    }

      if(length(featureSigID(MAIT.object))==0){
          stop("No significant features found in the MAIT object. Make sure that functions peakAnnotation and spectralSignFeatures were launched")
        }

    if(length(MAIT.object@FeatureData@masses)==0 & length(rawData(MAIT.object))==0){
      stop("No peak masses found in the MAIT object")
    }

    
      if (is.null(database)) {
                        identMetEnv<-new.env()
                data(MAITtables,envir=identMetEnv)
                Database<-get("Database",envir=identMetEnv)

        dataBase <- Database
       	cat("WARNING: No input database table was given. Selecting default MAIT database...",fill=TRUE)


      }else{
        
        dataBase<-read.csv(paste(database,".csv",sep=""),sep=",",header=TRUE)
	dataBase <- as.matrix(dataBase[order(dataBase[,1]),])
      }
        parameters <- list(peakTolerance,
                           database,
                           polarity)
        names(parameters) <- c("peakTolerance",
                               "database",
                               "polarity")

       MAIT.object@RawData@parameters@identifyMetabolites <- parameters
       writeParameterTable(parameters(MAIT.object),folder=resultsPath(MAIT.object))
    
##########################################################################################################################################################
#
#   Getting the necessary data from the MAIT object.
#
##########################################################################################################################################################

    
	signSpectra <-featureSigID(MAIT.object)
	sigPeaksTable <- sigPeaksTable(MAIT.object,printCSVfile=FALSE)
	resultsPath <- resultsPath(MAIT.object)

##########################################################################################################################################################
#
#   Defining the objects to save the results and the factors to perform the statistical tests.
#
##########################################################################################################################################################

    

    
    if(length(rawData(MAIT.object))==0){
      
      	Search <- matrix(nrow=1,ncol=8)
	colnames(Search) <- c("Query Mass","Database Mass (neutral mass)","rt","Adduct","Name","spectra","Biofluid","ENTRY")
        
      }else{
        
        Search <- matrix(nrow=1,ncol=9)
	colnames(Search) <- c("Query Mass","Database Mass (neutral mass)","rt","Isotope","Adduct","Name","spectra","Biofluid","ENTRY")
        
       }
    
    	H.mass <- 1.00794
    
#	peakList <- getPeaklist(MAIT.object)	
	

	temp <- getScoresTable(MAIT.object,getExtendedTable=TRUE)
        peakList <- temp$extendedTable
        spec <- temp$spectraID
        
    
        #signPeaklist <- peakList[as.numeric(peakList$pcgroup)%in%signSpectra,]
        signPeaklist <- peakList[spec%in%signSpectra,] 
    
	aux <- signPeaklist[c(-grep("[M+1]",signPeaklist$isotope,fixed=TRUE),-grep("[M+2]",signPeaklist$isotope,fixed=TRUE)),]
	aux <- aux[-which(aux$isotope==""),]

    if(length(rawData(MAIT.object))==0){

      	peaksIsTP <- matrix(nrow=1,ncol=dim(aux)[2]+2*length(classes(MAIT.object)))
	colnames(peaksIsTP) <- c(colnames(aux)[3:(dim(aux)[2]-1)],"p.adj","p",colnames(sigPeaksTable)[dim(sigPeaksTable)[2]-1-length(classes(MAIT.object))*2],colnames(sigPeaksTable)[(dim(sigPeaksTable)[2]-2*length(classes(MAIT.object))+1):dim(sigPeaksTable)[2]])

    }else{
    
	peaksIsTP <- matrix(nrow=1,ncol=dim(aux)[2]-7+2*length(classes(MAIT.object)))
	colnames(peaksIsTP) <- c(colnames(aux)[8:(dim(aux)[2]-3)],"p.adj","p",colnames(sigPeaksTable)[dim(sigPeaksTable)[2]-length(classes(MAIT.object))*2],colnames(sigPeaksTable)[(dim(sigPeaksTable)[2]-2*length(classes(MAIT.object))+1):dim(sigPeaksTable)[2]])
      }
  
  
	auxAd <- signPeaklist[which(signPeaklist$adduct!=""),]
	auxAd <- auxAd[order(as.numeric(auxAd$pcgroup)),]
	

	cat("Metabolite identification initiated",fill=TRUE)
        cat('\n % Metabolite identification in progress: ')
        lp <- -1
    
##########################################################################################################################################################
#
#   Each one of the significant features are looked up into the database.
#
##########################################################################################################################################################
    

	for (i in (1:dim(sigPeaksTable)[1])){
		spectra<-sigPeaksTable$pcgroup[i]
		index<-i

##########################################################################################################################################################
#
#   The search is adapted depending on the polarity of the raw data given by the input parameter polarity.
#
##########################################################################################################################################################

                
		if(polarity=="positive"){
			neutralPeak <- sigPeaksTable$mz[i]-H.mass
		}
                
		if(polarity=="negative"){
			neutralPeak <- sigPeaksTable$mz[i]+H.mass
		}

                
##########################################################################################################################################################
#
#   Funtion SearchCand looks up for a single peak into the database and returns all the suitable candidates.
#
##########################################################################################################################################################

		
		SearchPeak <- SearchCand(candidate=neutralPeak,dataBase=dataBase,peakTolerance=peakTolerance)

##########################################################################################################################################################
#
#   Once the search is performed, the results are attached to the identification table
#
##########################################################################################################################################################

                
		if(SearchPeak$unk==0){
			if (length(SearchPeak$SearchCand)==5){
				ref <- 1
                                
				SearchPeak$SearchCand <- matrix(SearchPeak$SearchCand,nrow=1)
			}else{
				ref <- dim(SearchPeak$SearchCand)[1]
			}
			for (k in c(1:ref)){

                          if(length(rawData(MAIT.object))==0){
                          
				Search <- rbind(Search,cbind(round(sigPeaksTable$mz[index],5),as.numeric(SearchPeak$SearchCand[k,4]),round(sigPeaksTable$rt[index],2),sigPeaksTable$adduct[index],as.character(SearchPeak$SearchCand[k,2]),spectra,as.character(SearchPeak$SearchCand[k,5]),as.character(SearchPeak$SearchCand[k,1])))
                                
                               if(length(classes(MAIT.object))>=2){
                                
			#	added <- cbind(round(sigPeaksTable[index,3:(dim(sigPeaksTable)[2]-5-length(classes(MAIT.object))*2)],0),sigPeaksTable$P.adjust[index],sigPeaksTable$p[index],sigPeaksTable[index,dim(sigPeaksTable)[2]-length(classes(MAIT.object))*2],sigPeaksTable[index,(dim(sigPeaksTable)[2]-2*length(classes(MAIT.object))+1):dim(sigPeaksTable)[2]])
				added <- cbind(sigPeaksTable[index,3:(dim(sigPeaksTable)[2]-5-length(classes(MAIT.object))*2)],sigPeaksTable$P.adjust[index],sigPeaksTable$p[index],sigPeaksTable[index,grep("Fisher.",colnames(sigPeaksTable))],sigPeaksTable[index,(dim(sigPeaksTable)[2]-2*length(classes(MAIT.object))+1):dim(sigPeaksTable)[2]])

                                 

                              }else{
                                
                                added <- cbind(sigPeaksTable[index,3:(dim(sigPeaksTable)[2]-5-length(classes(MAIT.object))*2)],sigPeaksTable$P.adjust[index],sigPeaksTable$p[index],sigPeaksTable$Fisher.Test[index],sigPeaksTable[index,(dim(sigPeaksTable)[2]-2*length(classes(MAIT.object))+1):dim(sigPeaksTable)[2]])
                                
                              }
				colnames(added) <- colnames(peaksIsTP)
				peaksIsTP <- rbind(peaksIsTP,added)

                           }else{

                                				Search <- rbind(Search,cbind(round(sigPeaksTable$mz[index],5),as.numeric(SearchPeak$SearchCand[k,4]),round(sigPeaksTable$rt[index],2),sigPeaksTable$isotopes[index],sigPeaksTable$adduct[index],as.character(SearchPeak$SearchCand[k,2]),spectra,as.character(SearchPeak$SearchCand[k,5]),as.character(SearchPeak$SearchCand[k,1])))
				#added <- cbind(round(sigPeaksTable[index,8:(dim(sigPeaksTable)[2]-6-length(classes(MAIT.object))*2)],0),sigPeaksTable$P.adjust[index],sigPeaksTable$p[index],sigPeaksTable$Fisher.Test[index],sigPeaksTable[index,(dim(sigPeaksTable)[2]-2*length(classes(MAIT.object))+1):dim(sigPeaksTable)[2]])

                                added <- cbind(sigPeaksTable[index,8:(dim(sigPeaksTable)[2]-6-length(classes(MAIT.object))*2)],sigPeaksTable$P.adjust[index],sigPeaksTable$p[index],sigPeaksTable[index,dim(sigPeaksTable)[2]-length(classes(MAIT.object))*2],sigPeaksTable[index,(dim(sigPeaksTable)[2]-2*length(classes(MAIT.object))+1):dim(sigPeaksTable)[2]])
				colnames(added) <- colnames(peaksIsTP)
				peaksIsTP <- rbind(peaksIsTP,added)
                      }
			}
                        
		}else{
			ref <- 1
                     if(length(rawData(MAIT.object))==0){

               			Search <- rbind(Search,cbind(round(sigPeaksTable$mz[index],5),"Unknown",round(sigPeaksTable$rt[index],2),sigPeaksTable$adduct[index],"Unknown",spectra,as.character(SearchPeak$SearchCand[5]),as.character(SearchPeak$SearchCand[1])))
#			added <- cbind(round(sigPeaksTable[index,3:(dim(sigPeaksTable)[2]-5-length(classes(MAIT.object))*2)],0),sigPeaksTable$P.adjust[index],sigPeaksTable$p[index],sigPeaksTable[index,dim(sigPeaksTable)[2]-1-length(classes(MAIT.object))*2],sigPeaksTable[index,(dim(sigPeaksTable)[2]-2*length(classes(MAIT.object))+1):dim(sigPeaksTable)[2]])
                       	added <- cbind(sigPeaksTable[index,3:(dim(sigPeaksTable)[2]-5-length(classes(MAIT.object))*2)],sigPeaksTable$P.adjust[index],sigPeaksTable$p[index],sigPeaksTable[index,grep("Fisher.",colnames(sigPeaksTable))],sigPeaksTable[index,(dim(sigPeaksTable)[2]-2*length(classes(MAIT.object))+1):dim(sigPeaksTable)[2]])
			colnames(added) <- colnames(peaksIsTP)
			peaksIsTP <- rbind(peaksIsTP,added)


                     }else{
                       
			Search <- rbind(Search,cbind(round(sigPeaksTable$mz[index],5),"Unknown",round(sigPeaksTable$rt[index],2),sigPeaksTable$isotopes[index],sigPeaksTable$adduct[index],"Unknown",spectra,as.character(SearchPeak$SearchCand[5]),as.character(SearchPeak$SearchCand[1])))

if(length(classes(MAIT.object))>=2){
                        
			added <- cbind(sigPeaksTable[index,8:(dim(sigPeaksTable)[2]-6-length(classes(MAIT.object))*2)],sigPeaksTable$P.adjust[index],sigPeaksTable$p[index],sigPeaksTable[index,grep("Fisher.",colnames(sigPeaksTable))],sigPeaksTable[index,(dim(sigPeaksTable)[2]-2*length(classes(MAIT.object))+1):dim(sigPeaksTable)[2]])

                      }else{
                        added <- cbind(sigPeaksTable[index,8:(dim(sigPeaksTable)[2]-6-length(classes(MAIT.object))*2)],sigPeaksTable$P.adjust[index],sigPeaksTable$p[index],sigPeaksTable$Fisher.Test[index],sigPeaksTable[index,(dim(sigPeaksTable)[2]-2*length(classes(MAIT.object))+1):dim(sigPeaksTable)[2]])
                      }
			colnames(added) <- colnames(peaksIsTP)
			peaksIsTP <- rbind(peaksIsTP,added)

                      }
		}
                           perc <-round(i/dim(sigPeaksTable)[1]*100)
           
           if ((perc %% 10 == 0) && (perc != lp)) { cat(perc,' '); lp <- perc }

              }
	Search <- Search[-1,]
	peaksIsTP <- peaksIsTP[-1,]
	peaksIsTP <- as.matrix(peaksIsTP)
        cat('\n Metabolite identification finished')
        row.names(peaksIsTP) <- 1:dim(peaksIsTP)[1]
      	row.names(Search) <- 1:dim(Search)[1]
        peaksIsTP.df <- as.data.frame(peaksIsTP)
        Search.df <- as.data.frame(Search)


    
##########################################################################################################################################################
#
#   The identification table is done and it is saved in the provided MAIT object. If the parameter printCSVfile is set to TRUE, the table is also printed
#   in the directory (working directory)/Tables/metaboliteTable.csv
#
##########################################################################################################################################################

    
    
	ID <- c(1:dim(Search)[1])
	metaboliteTable <- matrix(nrow=dim(Search)[1])
	rownames(metaboliteTable) <- rep("",dim(metaboliteTable)[1])
	rownames(Search) <- rep("",dim(metaboliteTable)[1])
	rownames(peaksIsTP) <- rep("",dim(peaksIsTP)[1])
        metaboliteTable <- merge(Search.df,peaksIsTP.df,by="row.names")[,-1]
        metaboliteTable[,1] <- round(as.numeric(as.character(metaboliteTable[,1])),4)
        #aux <- metaboliteTable
        #aux<-apply(aux,2,unlist)
        #aux[,-c(which(colnames(metaboliteTable)%in%c("Isotope","Adduct","Name","Biofluid","ENTRY")),2,dim(metaboliteTable)[2])] <- apply(aux[,-c(which(colnames(metaboliteTable)%in%c("Isotope","Adduct","Name","Biofluid","ENTRY")),2,dim(metaboliteTable)[2])],2,as.numeric)

    
##########################################################################################################################################################
#
#   Reordering the output table.
#
##########################################################################################################################################################

if(length(rawData(MAIT.object))==0){
  
ind1 <- 9:(dim(metaboliteTable)[2]-2*length(classes(MAIT.object))-3)
ind2<-(dim(metaboliteTable)[2]-2*length(classes(MAIT.object))-2):dim(metaboliteTable)[2]
tem<-metaboliteTable[,ind2]
metaboliteTable[,(dim(metaboliteTable)[2]+1-length(ind1)):dim(metaboliteTable)[2]] <- metaboliteTable[,ind1]
colnames(metaboliteTable)[(dim(metaboliteTable)[2]+1-length(ind1)):dim(metaboliteTable)[2]] <- colnames(metaboliteTable)[ind1]
metaboliteTable[,(dim(metaboliteTable)[2]+1-length(ind1)-length(ind2)):(dim(metaboliteTable)[2]-length(ind1))] <- tem
colnames(metaboliteTable)[(dim(metaboliteTable)[2]+1-length(ind1)-length(ind2)):(dim(metaboliteTable)[2]-length(ind1))] <- colnames(tem)

}else{

  ind1 <- 10:(dim(metaboliteTable)[2]-2*length(classes(MAIT.object))-3)
ind2<-(dim(metaboliteTable)[2]-2*length(classes(MAIT.object))-2):dim(metaboliteTable)[2]
tem<-metaboliteTable[,ind2]
metaboliteTable[,(dim(metaboliteTable)[2]+1-length(ind1)):dim(metaboliteTable)[2]] <- metaboliteTable[,ind1]
colnames(metaboliteTable)[(dim(metaboliteTable)[2]+1-length(ind1)):dim(metaboliteTable)[2]] <- colnames(metaboliteTable)[ind1]
metaboliteTable[,(dim(metaboliteTable)[2]+1-length(ind1)-length(ind2)):(dim(metaboliteTable)[2]-length(ind1))] <- tem
colnames(metaboliteTable)[(dim(metaboliteTable)[2]+1-length(ind1)-length(ind2)):(dim(metaboliteTable)[2]-length(ind1))] <- colnames(tem)

        #temp<-aux[,c((dim(aux)[2]-2):dim(aux)[2])]
       # spec <- as.numeric(as.character(aux[,7]))
        #biof <- aux[,8]
        #entry <- aux[,9]
        #temp2 <- aux[,10:(dim(aux)[2]-3)]
        #aux2<-aux    
        #aux2[,13:dim(aux2)[2]] <- temp2
        #colnames(aux2)[13:dim(aux2)[2]] <- colnames(temp2)
        #aux2[,10:12] <- temp
        #colnames(aux2)[10:12] <- colnames(temp)

#        metaboliteTable <- as.data.frame(aux2)
      }
    
        if(printCSVfile==TRUE){    
                                    	
	   if(!file.exists(paste(resultsPath,"Tables",sep="/"))){
		
		   dir.create(paste(resultsPath,"Tables",sep="/"))
                   
                   
                   
                    write.table(metaboliteTable,paste(paste(resultsPath(MAIT.object),"Tables","metaboliteTable",sep="/"),".csv",sep=""),col.names=NA,row.names=TRUE,sep=",")
		
	   }else{ 
		
		   cat(" " ,fill=TRUE)
		#   cat(paste("Warning: Folder",paste(resultsPath,"Tables",sep="/"),"already exists. Possible file overwritting.",sep=" "),fill=TRUE)
           	   warning(paste("Folder",paste(resultsPath,"Tables",sep="/"),"already exists. Possible file overwritting.",sep=" "))
                   write.table(metaboliteTable,paste(paste(resultsPath(MAIT.object),"Tables","metaboliteTable",sep="/"),".csv",sep=""),col.names=NA,row.names=TRUE,sep=",")


		
	   }
           
        }

        MAIT.object@FeatureInfo@metaboliteTable <- metaboliteTable

    
	return(MAIT.object)
    
}
