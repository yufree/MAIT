setMethod(f="loadings",signature="MAIT",function(object,type){
  if(!(exists("type"))){type<-"none"}
  
  out <- switch(type,
         "none" = models(object),
         "PCA" = pcaLoadings(object),
         "PLS" = plsLoadings(object))
  return(out)
})
