setMethod(f="scores",signature="MAIT",function(object,type="none"){
if(!(exists("type"))){type<-"none"}

 out <- switch(type,
         "none" = object@FeatureData@scores,
         "PCA" = pcaScores(object),
         "PLS" = plsScores(object))
  return(out)
})
