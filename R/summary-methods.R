setMethod(f="summary",signature="MAIT",function(object){
  show(object)  
  if(!(dim(ovClassifRatioTable(object))[1]==1)){            
                  cat(paste("The Classification using",parameters(object)@classification[[2]],"training samples and",parameters(object)@classification[[1]],"Iterations gave the results:","\n",sep=" "),fill=TRUE)
  show(ovClassifRatioTable(object));
  }
  cat("\n",fill=TRUE)
   cat("Parameters of the analysis:",fill=TRUE)
writeParameterTable(parameters(object),folder=resultsPath(object))}
)
