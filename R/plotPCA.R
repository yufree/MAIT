########################################################################################################################################################## 
##########################################################################################################################################################
##                                                      |
##   Metabolite Automatic Identification Toolkit (MAIT) |
##                                                      |
##                                                      |
##   written by Francesc Fernández Albert               |
##   contact mail: francesc.fernandez.albert@upc.edu    |
##   date: 7/29/2013                                    |
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



plotPCA<-function (MAIT.object = NULL, Log = FALSE, center = TRUE, scale = TRUE, plot3d=TRUE) 
{

    if (is.null(MAIT.object)) {
        stop("No input MAIT object file was given")
    }
    if (length(featureSigID(MAIT.object)) == 0) {
        stop("No significant features found in the MAIT object. Make sure that functions peakAnnotation and peakAggregation were launched")
    }
    
    parameters <- list(Log, center, scale)
    names(parameters) <- c("PCA data logarithm", "PCA data centered", 
        "PCA data scaled")
    MAIT.object@RawData@parameters@plotPCA <- parameters
    writeParameterTable(parameters(MAIT.object), folder = resultsPath(MAIT.object))
    data <- scores(MAIT.object)
    clases <- classes(MAIT.object)
    classNum <- classNum(MAIT.object)
    xsaFA <- MAIT.object@RawData@data
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
    if (!file.exists(paste(resultsPath, "PCA_Scoreplots", sep = "/"))) {
        dir.create(paste(resultsPath, "PCA_Scoreplots", sep = "/"))
    }
    else {
        cat(" ", fill = TRUE)
        cat(paste("Warning: Folder", paste(resultsPath, "PCA_Scoreplots", 
            sep = "/"), "already exists. Possible file overwritting.", 
            sep = " "), fill = TRUE)
    }

    model <- prcomp(data)

    
    png(paste(resultsPath, "PCA_Scoreplots/Scoreplot_PC12.png", 
        sep = "/"))   
    
    plot(x = model$x[, 1], y = model$x[, 2], col = cols,pch=16,cex=0.8, xlab=paste("PC1 (",round(summary(model)$importance[2,1],3),")",sep=""),ylab=paste("PC2 (",round(summary(model)$importance[2,2],3),")",sep=""))
    legend("topleft", legend = clases, text.col = textCols)
    dev.off()
    png(paste(resultsPath, "PCA_Scoreplots/Scoreplot_PC23.png", 
        sep = "/"))
    plot(x = model$x[, 2], y = model$x[, 3], col = cols,pch=16,cex=0.8, xlab=paste("PC2 (",round(summary(model)$importance[2,2],3),")",sep=""),ylab=paste("PC3 (",round(summary(model)$importance[2,3],3),")",sep=""))
    legend("topleft", legend = clases, text.col = textCols)
    dev.off()
    png(paste(resultsPath, "PCA_Scoreplots/Scoreplot_PC13.png", 
        sep = "/"))
    plot(x = model$x[, 1], y = model$x[, 3], col = cols,pch=16,cex=0.8, xlab=paste("PC1 (",round(summary(model)$importance[2,1],3),")",sep=""),ylab=paste("PC3 (",round(summary(model)$importance[2,3],3),")",sep=""))
    legend("topleft", legend = clases, text.col = textCols)
    dev.off()
    PC1 <- model$x[, 1]
    PC2 <- model$x[, 2]
    PC3 <- model$x[, 3]

    if  (plot3d == TRUE){
    require(rgl)
    PCAplot3d(x = PC1, y = PC2, z = PC3, cols = cols)
  }
    MAIT.object@FeatureData@pcaModel <- list(model)

    return(MAIT.object)
    
}


PCAplot3d<-function(z, x, y, cols,axes=TRUE,new=TRUE)

{xr<-range(x) 
x01<-(x-xr[1])/(xr[2]-xr[1]) 
yr<-range(y) 
y01<-(y-yr[1])/(yr[2]-yr[1]) 
zr<-range(z) 
z01<-(z-zr[1])/(zr[2]-zr[1])

if(new) rgl::rgl.clear() 
if(axes) 
        {xlab<-pretty(x) 
        ylab<-pretty(y) 
        zlab<-pretty(z) 
        xat<-(xlab-xr[1])/(xr[2]-xr[1]) 
        yat<-(ylab-yr[1])/(yr[2]-yr[1]) 
        zat<-(zlab-zr[1])/(zr[2]-zr[1]) 
        rgl::rgl.lines(c(0,1.1),0,0) 
        rgl::rgl.lines(0,c(0,1.1),0) 
        rgl::rgl.lines(0,0,c(0,1.1))	
        rgl::rgl.texts(xat,-.05,-.05,xlab) 
        rgl::rgl.texts(-.05,yat,-.05,ylab) 
        rgl::rgl.texts(-.05,-.05,zat,zlab) 
        rgl::rgl.texts(c(0.5,-.15,-.15),c(-.15,.5,-.15),c(-.15,-.15,.5), 
                c(deparse(substitute(x)),deparse(substitute(y)),deparse(substitute(z)))) 
        } 

rgl::rgl.spheres(x01,y01,z01,.01,color=(cols)) 
} 
