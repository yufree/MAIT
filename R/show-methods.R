setMethod(f="show",signature="MAIT",
          function(object) {

            if(length(rawData(object))==0 & dim(object@FeatureData@scores)[1]==1 & length(object@FeatureData@masses)==0){
              
                cat("This MAIT object does not have any data information",fill=TRUE)

              }else{
                


      if (method(object)=="None"|method(object)==""){
  
                vars <- "peaks"
                if(dim(scores(object))[1]!=1){
                cat(paste("A MAIT object built of",sum(classNum(object)),"samples and",nrow(scores(object)),paste(vars,".",sep=""),"No peak aggregation technique has been applied"),fill=TRUE);
              }else{
                cat(paste("A MAIT object built of",sum(classNum(object)),"samples"),fill=TRUE);
              }
                
      if(!(length(featureSigID(object))==0)){                
               cat(paste(length(featureSigID(object)),"of these peaks are statistically significant",sep=" "),fill=TRUE);
             } 

                          for (i in c(1:length(classes(object)))){
                          cat(paste("The object contains",classNum(object)[i],"samples of class",classes(object)[i],sep=" "),fill=TRUE)
                        }

                
          }else{   
          vars <- "spectra"
         
          NcolData<-ncol(scores(object))
          NrowData<-nrow(scores(object))
 
        cat(paste("An MAIT object built of",NcolData,"samples and",NrowData,vars,sep=" "),fill=TRUE);
                if(!(length(featureSigID(object))==0)){
                
               cat(paste(length(featureSigID(object)),"of these spectra are statistically significant",sep=" "),fill=TRUE);
             }


        cat(paste("Peak aggregation method used:",method(object),"\n",sep=" "),fill=FALSE);


                          for (i in c(1:length(classes(object)))){
                          cat(paste("The object contains",classNum(object)[i],"samples of class",classes(object)[i],sep=" "),fill=TRUE)
                        }

	}
       }
          
          }
)
      




