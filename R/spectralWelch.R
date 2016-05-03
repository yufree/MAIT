########################################################################################################################################################## 
##########################################################################################################################################################
##                                                      |
##   Metabolite Automatic Identification Toolkit (MAIT) |
##                                                      |
##                                                      |
##   written by Francesc Fernández Albert               |
##   contact mail: francesc.fernandez.albert@upc.edu    |
##   date: 10/15/2012                                   |
##                                                      |
##   SISBIO Group. ESAII Department                     |
##   Universitat Politècnica de Catalunya               |
##   ___________________________________________________|
##
##    
##   Function spectralAnova takes a numerical matrix and performs an ANOVA test for each of its variables using the information of a MAIT object. P-value threshold could be
##   defined to build an output table made of significant features only. Fisher LSD analysis is also performed in order to know which are the class differences in the
##   statistically significant features.
##   
##   MAIT.object:  Object containing the spectral data after having applied the function peakAggregation over it. T-tests will be performed on each of its features.
##   pvalue:       P-values having a lower value than pvalue are significant and therefore are included in the output (set to 0.05 by default).
##   bonferroni:   If it is set to TRUE, Bonferroni multiple testing correction is performed. Pvalue threshold is taken after the correction
##   printCSVfile: If it is set to TRUE, a table containing the results of the statistical tests is printed under the name (working directory)/Tables/significativeFeatures.csv
##                 through the signPeaksTable function
##
##
## All code copyright (c) 2013 UPC/UB
## All accompanying written materials copyright (c) 2013 UPC/UB
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
########################################################################################################################################################## 
##########################################################################################################################################################


spectralWelch <- function(MAIT.object=NULL,
                             pvalue=0.05,
                             p.adj="none",
                             printCSVfile=TRUE){
  
    if (is.null(MAIT.object)) {
        stop("No input MAIT object was given")
    }


        parameters <- list(pvalue,
                           p.adj)
        names(parameters) <- c("Welch pvalue",
                           "Welch p.adj")


       MAIT.object@RawData@parameters@sigFeatures <- parameters
       writeParameterTable(parameters(MAIT.object),folder=resultsPath(MAIT.object))
    
##########################################################################################################################################################
#
#   Getting the necessary data from the MAIT object.
#
##########################################################################################################################################################

	
        data <- scores(MAIT.object)
        clases <- classes(MAIT.object)
        classNum <- classNum(MAIT.object)
        #xsaFA <- rawData(MAIT.object)
        resultsPath <- resultsPath(MAIT.object)
#	peakList <- getPeaklist(MAIT.object)

##########################################################################################################################################################
#
#   Creating the objects to save the results and the factors to perform the statistical tests.
#
##########################################################################################################################################################

    
#	peakList <- peakList[order(as.numeric(peakList$pcgroup)),]

        auxs<-getScoresTable(MAIT.object=MAIT.object,getExtendedTable=TRUE)
        peakList <- auxs$extendedTable
        spec <- auxs$spectraID
    
	classes <- c(rep(clases[1],classNum[1]),rep(clases[2],classNum[2]))
	numbers <- classNum
	Bonferroni <- matrix(ncol=as.numeric(peakList$pcgroup[dim(peakList)[1]]),nrow=1)
	names(numbers) <- clases
	TTs <- matrix(ncol=1,nrow=dim(data)[1])	
	colnames(TTs) <- paste(clases[1],"&",clases[2],sep="")
#	Tresults <- matrix(ncol=1,nrow=as.numeric(peakList$pcgroup[dim(peakList)[1]]))
    	Tresults <- matrix(ncol=1,nrow=spec[dim(peakList)[1]])
	group <- as.factor(c(rep(clases[1],numbers[1]),rep(clases[2],numbers[2])))

    
##########################################################################################################################################################
#
#   Loop over all the features in order to perform a statistical test over each variable.
#
##########################################################################################################################################################

    
    for (i in c(1:(dim(data)[1]))){
		numbers <- classNum
		names(numbers) <- clases
		lmdata <- as.vector(t(data[i,]))
		model <- lm(lmdata~group)
		if (length(which(diff(data)!=0))){
			ttest <- t.test(as.vector(t(data[i,]))~group,var.equal=FALSE)
			TTs[i] <- ttest$p.value
			Tresults[i] <- ttest$statistic
		}else{
			TTs[i] <- 1
		}
	}

    
##########################################################################################################################################################
#
#   The p-values of the statistical tests are saved in the MAIT object.
#
##########################################################################################################################################################

    
        MAIT.object@FeatureData@pvalues <- TTs
	 p.corr <- p.adjust(TTs,method=p.adj)


##########################################################################################################################################################
#
#   Depending on the value of the parameter bonferroni, the significant features are chosen by corrected p-values or not.
#
##########################################################################################################################################################

    
	if (p.adj!="none"){
		
		index <- which(p.corr<=pvalue)

		
	}else{
		
		index<-which(TTs<=pvalue)
		
	}
    
        MAIT.object@FeatureData@pvaluesCorrection <- p.adj
        MAIT.object@FeatureData@featureSigID<-index

    

   # if(printCSVfile==TRUE){
    
   #   aux <- signPeaksTable(MAIT.object,printCSVfile=printCSVfile)

    #}

    
  return(MAIT.object)  
}

    
