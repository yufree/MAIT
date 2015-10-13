########################################################################################################################################################## 
##########################################################################################################################################################
##                                                      |
##   Metabolite Automatic Identification Toolkit (MAIT) |
##                                                      |
##                                                      |
##   written by Francesc Fernández Albert               |
##   contact mail: francesc.fernandez.albert@upc.edu    |
##   date: 7/29/2013                                   |
##                                                      |
##   SISBIO Group. ESAII Department                     |
##   Technical University of Catalonia                  |
##   Nutrition Department                               |
##   University of Barcelona                            |
##   ___________________________________________________|
##
##    
##   signBoxplot function builds a png file containing a boxplot for each significative feature present in a MAIT object. These files are saved in a folder
##   called Boxplots inside the project directory.
##    
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


plotBoxplot <- function(MAIT.object=NULL){

    if (is.null(MAIT.object)) {
        stop("No input MAIT object file was given")
    }

    if(length(featureSigID(MAIT.object))==0){
          stop("No significant features found in the MAIT object. Make sure that functions peakAnnotation and peakAggregation were launched")
        }

    
        data <- scores(MAIT.object)
	index <- featureSigID(MAIT.object)
        class <- classes(MAIT.object)
        classNum <- classNum(MAIT.object)
        resultsPath <- resultsPath(MAIT.object)

  
	clases <- matrix(nrow=1)
	for (i in c(1:length(class))){
		clases <- c(clases,rep(class[i],classNum[i]))
	}
	clases <- clases[-1]
	clases <- as.factor(clases)
	aux <- t(data[index,])
    				if(!file.exists(paste(resultsPath,"Boxplots",sep="/"))){
					
					dir.create(paste(resultsPath,"Boxplots",sep="/"))

				}else{

					cat(" " ,fill=TRUE)
					cat(paste("Warning: Folder",paste(resultsPath,"Boxplots",sep="/"),"already exists. Possible file overwritting.",sep=" "),fill=TRUE)

				}

	for (i in c(1:length(index))){
		png(paste(paste(resultsPath,"Boxplots/Boxplot_spectra_",sep="/"),index[i],".png",sep=""))
		boxplot(aux[,i]~clases)
		title(paste("Spectra",index[i],sep=""))
		dev.off()
	}
}
