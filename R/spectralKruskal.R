########################################################################################################################################################## 
##########################################################################################################################################################
##                                                      |
##   Metabolite Automatic Identification Toolkit (MAIT) |
##                                                      |
##                                                      |
##   written by Francesc Fernández Albert               |
##   contact mail: francesc.fernandez.albert@upc.edu    |
##   date: 1/10/2014                                    |
##                                                      |
##   SISBIO Group. ESAII Department                     |
##   Technical University of Catalonia                  |
##   Nutrition Department                               |
##   University of Barcelona                            |
##                                                      |
##   SISBIO Group. ESAII Department                     |
##   Universitat Politècnica de Catalunya               |
##   ___________________________________________________|
##
##    
##   Function spectralAnova takes a numerical matrix and performs an ANOVA test for each of its variables using the information of a MAIT object. P-value threshold
##   could be defined to build an output table made of significant features only. Fisher LSD analysis is also performed in order to know which are the class
##   differences in the statistically significant features.
##   
##   MAIT.object:  Object containing the spectral data after having applied the function peakAggregation over it. Non-parametric Krustal-Wallis tests will be performed on each of its
##                 features.
##   pvalue:       P-values having a lower value than pvalue are significant and therefore are included in the output (set to 0.05 by default).
##   bonferroni:   If it is set to TRUE, Bonferroni multiple testing correction is performed. Pvalue threshold is taken after the correction.
##   printCSVfile: If it is set to TRUE, a table containing the results of the statistical tests is printed under the name
##                 (working directory)/Tables/significativeFeatures.csv is printed through the signPeaksTable function
##                 
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


spectralKruskal <- function (pvalue=0.05,
                           p.adj="none",
                           MAIT.object=NULL,
                           printCSVfile=TRUE)
{

    if (is.null(MAIT.object)) {
        stop("No input MAIT object was given")
    }

        parameters <- list(pvalue,
                            p.adj)
        names(parameters) <- c("Kruskal-Wallis p-value",
                           "Kruskal-Wallis  p.adj")

  
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
        resultsPath <- resultsPath(MAIT.object)
        auxs<-getScoresTable(MAIT.object=MAIT.object,getExtendedTable=TRUE)
        peakList <- auxs$extendedTable
        
        
##########################################################################################################################################################
##       
##   Variable Fgroups contains the factor groups needed to calculate the test.
##
########################################################################################################################################################## 

        

	Fgroups <- matrix(nrow=1)
        Fgroups <- rep(clases,classNum)    
	Fgroups <- as.factor(Fgroups[order(Fgroups)])
        
##########################################################################################################################################################
##       
##  TTs variable contains the p-values of all the tests
##
########################################################################################################################################################## 

        
        TTs <- matrix(ncol=1,nrow=dim(data)[1])
        Fisher <- matrix(ncol=1,nrow=dim(data)[1])
        
	Tresults <- matrix(ncol=1,nrow=as.numeric(peakList$pcgroup[dim(peakList)[1]]))
    
    
	for (i in c(1:dim(data)[1])){
        
	   lmdata <- as.vector(data[i,])
	   numbers <- classNum
	   names(numbers) <- clases
#	   model <- lm(lmdata~Fgroups)
	   an <- kruskal.test(x=lmdata,g=Fgroups)
	   TTs[i] <- an$p.value
        
	}

    MAIT.object@FeatureData@pvalues <- TTs
#    MAIT.object@FeatureData@LSDResults <- Fisher
    
##########################################################################################################################################################
##       
##  Only those p-values having values under pvalue are considered statistically relevant. 
##
##########################################################################################################################################################
        
	 p.corr <- p.adjust(TTs,method=p.adj)

##########################################################################################################################################################
##       
##  Variable index is a vector containing the variable's ID of the statistically significant variables of input data.
##
##########################################################################################################################################################

        
	if (p.adj!="none"){
		
		index <- which(p.corr<=pvalue)

		
	}else{
		
		index<-which(TTs<=pvalue)
		
	}
    
        MAIT.object@FeatureData@pvaluesCorrection <- p.adj
        MAIT.object@FeatureData@featureSigID<-index

    

  return(MAIT.object)  
}
