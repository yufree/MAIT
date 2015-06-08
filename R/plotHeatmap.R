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
##   signHeatmap function creates few heatmaps as png files considering different distances and p-values when it is applied on a MAIT object.
##
##   MAIT.object:  MAIT object containing the data after having gathered the significant features.
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


plotHeatmap<-function (MAIT.object = NULL) 
{
    if (is.null(MAIT.object)) {
        stop("No input MAIT object file was given")
    }
    if (length(featureSigID(MAIT.object)) == 0) {
        stop("No significant features found in the MAIT object. Make sure that functions peakAnnotation and peakAggregation were launched")
    }

    writeParameterTable(parameters(MAIT.object), folder = resultsPath(MAIT.object))
    data <- scores(MAIT.object)
    index <- featureSigID(MAIT.object)
    clases <- classes(MAIT.object)
    classNum <- classNum(MAIT.object)
    
    if(length(rawData(MAIT.object))!=0){
    xsaFA <- MAIT.object@RawData@data[[1]]
    names <- sampnames(xsaFA@xcmsSet)
    }else{
    names <- colnames(data)
    }
    
    resultsPath <- resultsPath(MAIT.object)
    names <- vector(length = sum(classNum))
    pvalues <- pvalues(MAIT.object)
    if(pvaluesCorrection(MAIT.object)=="Bonferroni"){
      pvalues <- p.adjust(pvalues)
    }
    distCor <- function(x) {
        distCor <- as.dist(1 - cor(t(x)))
    }
    hclustWard <- function(distCor) {
        hclustWard <- hclust(distCor, "ward")
    }
    rownames(data) <- as.character(c(1:dim(data)[1]))
    colnames(data) <- names
    cols <- matrix(nrow = 1)
    for (i in c(1:length(classNum))) {
        cols <- c(cols, rep(i, classNum[i]))
    }
    label <- 1:length(classNum)
    cols <- as.character(cols[-1])
    paletteCols <- colorRampPalette(c("cyan", "yellow"))(28)
    pvals <- c(0.05, 0.01, 0.001, 1e-04, 1e-05)

    cexRows <-   c(tanh(300/length(which(pvalues(MAIT.object)<=pvals[1]))),tanh(300/length(which(pvalues(MAIT.object)<=pvals[2]))),tanh(300/length(which(pvalues(MAIT.object)<=pvals[3]))),tanh(300/length(which(pvalues(MAIT.object)<=pvals[4]))),tanh(300/length(which(pvalues(MAIT.object)<=pvals[5]))))
    tanh(500/10000)
    points <- 0.2*c(length(which(pvalues(MAIT.object)<=pvals[1])),length(which(pvalues(MAIT.object)<=pvals[2])),length(which(pvalues(MAIT.object)<=pvals[3])),length(which(pvalues(MAIT.object)<=pvals[4])),length(which(pvalues(MAIT.object)<=pvals[5])))

    
    heights <- c(5020, 4096, 3072, 2048, 1024)
    if (!file.exists(paste(resultsPath, "Heatmaps", sep = "/"))) {
        dir.create(paste(resultsPath, "Heatmaps", sep = "/"))
    }else {
        cat(" ", fill = TRUE)
        cat(paste("Warning: Folder", paste(resultsPath, "Heatmaps", 
            sep = "/"), "already exists. Possible file overwritting.", 
            sep = " "), fill = TRUE)
    }
    for (i in c(1:length(pvals))) {
        if (length(which(pvalues <= pvals[i])) > 2) {
            png(paste(paste(resultsPath, "Heatmaps/Correlation_Distance_Heatmap_p", 
                sep = "/"), pvals[i], ".png", sep = ""), height = heights[i], 
                width = 1024, units = "px", pointsize = 40)
            if (nchar(names)[1] > 8) {
                colLength <- 0.45
            }else {
                colLength <- 0.8
            }

            heatmap.2(t(scale(t(data[which(pvalues <= pvals[i]), 
                ]), center = TRUE, scale = TRUE)),
                hclustfun = hclustWard, ColSideColors = cols, 
                trace = "none", distfun = distCor, cexCol = colLength, 
                cexRow = as.numeric(paste("0.", 8 - i, sep = "")),breaks=c(seq(-5,-1.5,by=1),seq(-1,1,by=0.1),seq(2,5,by=1)),col=paletteCols,symbreaks=TRUE)
            legend("topright", legend=clases, text.col = label, cex = as.numeric(paste("0.", 
                9 - i, sep = "")))
            dev.off()
            png(paste(paste(resultsPath, "Heatmaps/Euclidean_Distance_Heatmap_p", 
                sep = "/"), pvals[i], ".png", sep = ""), height = heights[i], 
                width = 1024, units = "px", pointsize = 40)
            heatmap.2(t(scale(t(data[which(pvalues <= pvals[i]), 
                ]), center = TRUE, scale = TRUE)),
                hclustfun = hclustWard, ColSideColors = cols, 
                trace = "none", cexCol = colLength, cexRow = as.numeric(paste("0.", 8 - i, sep = "")),breaks=c(seq(-5,-1.5,by=1)
                ,seq(-1,1,by=0.1),seq(2,5,by=1)),col=paletteCols,symbreaks=TRUE)
            legend("topright", legend=clases, text.col = label, cex = as.numeric(paste("0.", 
                9 - i, sep = "")))
            dev.off()
        }
        else {
            cat(paste("WARNING: There are not enough significant spectra for pvalue", 
                pvals[i], sep = " "), fill = TRUE)
        }
    }
}
