########################################################################################################################################################## 
##########################################################################################################################################################
##                                                      |
##   Metabolite Automatic Identification Toolkit (MAIT) |
##                                                      |
##                                                      |
##   written by Francesc Fernández Albert               |
##   contact mail: francesc.fernandez.albert@upc.edu    |
##   date: 1/09/2014                                    |
##                                                      |
##   SISBIO Group. ESAII Department                     |
##   Technical University of Catalonia                  |
##   Nutrition Department                               |
##   University of Barcelona                            |
##   ___________________________________________________|
##
##    
##   plotPCA function builds three png files containing the projection of the data over the first three principal components. Additionaly an interactive 3D
##   PCA plot is generated using the rgl package.
##
##   MAIT.object:     MAIT Object containing the data after having gathered the significant features.
##   Log:             If it is set to TRUE, the logarithm of the data is taken.
##   center:          If it is set to TRUE, the data is centered around its mean.
##   scale:           If it is set to TRUE, the data is scale to have unit variance.
##
## All code copyright (c) 2013 Francesc Fernandez
## All accompanying written materials copyright (c) 2013 Francesc Fernandez
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



plotPLS<-function (MAIT.object = NULL, Log = FALSE, center = TRUE, scale = TRUE, plot3d = TRUE) 
{

    if (is.null(MAIT.object)) {
        stop("No input MAIT object file was given")
    }
    if (length(featureSigID(MAIT.object)) == 0) {
        stop("No significant features found in the MAIT object. Make sure that functions peakAnnotation and peakAggregation were launched")
    }
    
    parameters <- list(Log, center, scale)
    names(parameters) <- c("PLS data logarithm", "PLS data centered", 
        "PLS data scaled")
    MAIT.object@RawData@parameters@plotPLS <- parameters
    writeParameterTable(parameters(MAIT.object), folder = resultsPath(MAIT.object))
    data <- scores(MAIT.object)
    clases <- classes(MAIT.object)
    classNum <- classNum(MAIT.object)
    resultsPath <- resultsPath(MAIT.object)
    index <- featureSigID(MAIT.object)
    cols <- matrix(nrow = 1)
    textCols <- matrix(nrow = 1)
    for (i in 1:length(clases)) {
        cols <- c(cols, rep(i, classNum[i]))
    }
    cols <- as.character(cols[-1])
    textCols <- 1:length(clases)
    if (Log == FALSE) {
        data <- (scale(t(data[index, ]), center = center, scale = scale))
    }
    else {
        data <- (scale(t(log10(data[index, ]+1)), center = center, 
            scale = scale))
    }
    if (!file.exists(paste(resultsPath, "PLS_Scoreplots", sep = "/"))) {
        dir.create(paste(resultsPath, "PLS_Scoreplots", sep = "/"))
    }
    else {
        cat(" ", fill = TRUE)
        cat(paste("Warning: Folder", paste(resultsPath, "PLS_Scoreplots", 
            sep = "/"), "already exists. Possible file overwritting.", 
            sep = " "), fill = TRUE)
    }
    png(paste(resultsPath, "PLS_Scoreplots/Scoreplot_PLS12.png", 
        sep = "/"))

    
    model <- train(x=data,y=as.factor(rep(clases,classNum)),method="pls")

if(model$finalModel$ncomp==1){
  warning("Just one component was found in the PLS")
      png(paste(resultsPath, "PLS_Scoreplots/Scoreplot_PLS1.png", sep = "/"))
      plot(y = matrix(model$finalModel$scores), x=c(1:length(matrix(model$finalModel$scores))),col = cols,pch=16,cex=0.8, ylab=paste("PC1 (",round((model$finalModel$Xvar/ model$finalModel$Xtotvar)[1],3),")",sep=""),xlab="Sample Index")
      legend("topright", legend = clases, text.col = textCols)
  dev.off()
}

    
if(model$finalModel$ncomp>1){
  
    png(paste(resultsPath, "PLS_Scoreplots/Scoreplot_PLS12.png", sep = "/"))
    
    plot(x = model$finalModel$scores[, 1], y = model$finalModel$scores[, 2], col = cols,pch=16,cex=0.8, xlab=paste("PC1 (",round((model$finalModel$Xvar/ model$finalModel$Xtotvar)[1],3),")",sep=""),ylab=paste("PC2 (",round((model$finalModel$Xvar/ model$finalModel$Xtotvar)[2],3),")",sep=""))
    legend("topleft", legend = clases, text.col = textCols)
    dev.off()
  }

    if(model$finalModel$ncomp>2){
      
       png(paste(resultsPath, "PLS_Scoreplots/Scoreplot_PLS23.png", sep = "/"))
 
    plot(x = model$finalModel$scores[, 2], y = model$finalModel$scores[, 3], col = cols,pch=16,cex=0.8, xlab=paste("PC2 (",round((model$finalModel$Xvar/ model$finalModel$Xtotvar)[2],3),")",sep=""),ylab=paste("PC3 (",round((model$finalModel$Xvar/ model$finalModel$Xtotvar)[3],3),")",sep=""))
    legend("topleft", legend = clases, text.col = textCols)
    dev.off()
    png(paste(resultsPath, "PLS_Scoreplots/Scoreplot_PLS13.png", sep = "/"))
    plot(x = model$finalModel$scores[, 1], y = model$finalModel$scores[, 3], col = cols,pch=16,cex=0.8, xlab=paste("PC1 (",round((model$finalModel$Xvar/ model$finalModel$Xtotvar)[1],3),")",sep=""),ylab=paste("PC3 (",round((model$finalModel$Xvar/ model$finalModel$Xtotvar)[3],3),")",sep=""))
    legend("topleft", legend = clases, text.col = textCols)
    dev.off()
    PC1 <- model$finalModel$scores[, 1]
    PC2 <- model$finalModel$scores[, 2]
    PC3 <- model$finalModel$scores[, 3]
  }

    if  (plot3d == TRUE){
    require(rgl)
    PCAplot3d(x = PC1, y = PC2, z = PC3, cols = cols)
  }

    
MAIT.object@FeatureData@plsModel <- list(model)
     
    return(MAIT.object)
    
}
