setMethod(f="model",signature="MAIT",function(x,type){if(is.null(type)==TRUE){
stop("No type was specified (PCA or PLS).")
}

  out <- switch(type,
         "PCA" = pcaModel(x),
         "PLS" = plsModel(x))

  return(out)})
