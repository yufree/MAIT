setMethod(f="getPeaklist",signature="MAIT",
          function(object) {
            if(class(rawData(object)[[1]])=="xsAnnotate"){
              return(CAMERA::getPeaklist(rawData(object)[[1]]))}else{
              }
            stop("Peak annotation has not been performed. Run function peakAnnotation first.")
          }
            )
