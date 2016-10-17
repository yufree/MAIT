########################################################################################################################################################## 
##########################################################################################################################################################
##                                                      |
##   Metabolite Automatic Identification Toolkit (MAIT) |
##                                                      |
##                                                      |
##   written by Francesc Fernández Albert               |
##   contact mail: francesc.fernandez.albert@upc.edu    |
##   date: 10/08/2013                                    |
##                                                      |
##   SISBIO Group. ESAII Department                     |
##   Technical University of Catalonia                  |
##   Nutrition Department                               |
##   University of Barcelona                            |
##   ___________________________________________________|
##
##
##   peakAnnotation function performs spectra building and peak annotation using the CAMERA package on a MAIT object, after applying the signalProcessing function.
##   The resultant xsAnnotate object is stored in a MAIT object. The input variables of the function are:
##
##   MAIT.object:       A MAIT object where function signalProcessing has already been applied. The output of the function is going to be an update of the same MAIT##                      object.
##   corr:              Correlation threshold value. If the correlation between two peaks is higher than this value, they are considered as peaks of the same       ##                      spectrum. Otherwise they are splitted in two different spectra.
##   perfwhm:           This parameter is used to group two peaks depending on their retention time. Two peaks are considered to be coeluted if their retention time##                      falls in a range defined as  Rt_med +/- FWHM * perfwhm. Where Rt_med is the retention time median whereas FWHM is the Full Width at         ##                      Half Maximum.
##                      Defined this way, perfwhm is the percentage of the width of the FWHM (Full Width at Half Maximum)
##   sigma:             Defining the coelution range as defined in the perfwhm variable, the FWHM is obtained by FWHM = SD * sigma, where SD is calculated          ##                      considering the peak as normally distributed. The multiplier of the standard deviation.
##   
##   adductTable:       Input table containing the adducts or fragments to be looked for in each spectrum.
##   printSpectraTable: If it is set to TRUE, a smaller version of the CAMERA getPeaklist table is printed.
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



peakAnnotation <- function(MAIT.object=NULL,
                           corrWithSamp=0.7,
                           perfwhm=0.6,
                           sigma=6,
                           adductTable=NULL,
                           printSpectraTable=TRUE,
                           corrBetSamp=0.75,
                           pval=0.05,
                           calcIso=TRUE,
                           calcCiS=TRUE,
                           calcCaS=TRUE,
                           graphMethod="hcs",
                           annotateAdducts=TRUE){


  
    if (is.null(MAIT.object)) {
        stop("No MAIT object was given")
    }

    if (length(MAIT.object@RawData@data)!=1) {
        stop("Raw data is not correct. Try running again signalProcessing function")
    }


        parameters <- list(corrWithSamp,
                          corrBetSamp,
                           perfwhm,
                           sigma,
                           adductTable,
                           pval,
                           calcIso,
                           calcCiS,
                           calcCaS,
                           graphMethod,
                           annotateAdducts)


     names(parameters) <- c("corrWithSamp",
                            "corrBetSamp",
                           "perfwhm",
                           "sigma",
                           "adductTable",
                            "peakAnnotation pvalue",
                            "calcIso",
                            "calcCiS",
                            "calcCaS",
                            "graphMethod",
                            "annotateAdducts")

  
       MAIT.object@RawData@parameters@peakAnnotation <- parameters
       writeParameterTable(parameters(MAIT.object),folder=resultsPath(MAIT.object))

  
##########################################################################################################################################################
#
#  If the adductTable input parameter is not defined, the function takes the default MAIT table for positive polarization as adductTable. If adductTable parameter is set to
#  "negAdducts", the default MAIT table for negative adducts is taken instead.
#
##########################################################################################################################################################


    if (is.null(adductTable)){
      
		cat("WARNING: No input adduct/fragment table was given. Selecting default MAIT table for positive polarity...",fill=TRUE)
                cat("Set adductTable equal to negAdducts to use the default MAIT table for negative polarity",fill=TRUE)
                peakAnnEnv<-new.env()
                data(MAITtables,envir=peakAnnEnv)
         adducts<-get(x="posAdducts",envir=peakAnnEnv)
       }else{

    if (adductTable=="negAdducts"){
    cat("adductTable has been set to negAdducts. The default MAIT adducts table for negative polarization is selected...",fill=TRUE)

                peakAnnEnv<-new.env()
                data(MAITtables,envir=peakAnnEnv)
      adducts<-get(x="negAdducts",envir=peakAnnEnv)
    }else{

         	adducts <- read.csv2(paste(adductTable,".csv",sep=""),dec=".",header=TRUE,sep=",")
       }
  }

    
##########################################################################################################################################################
#
#  Peak grouping (Spectra) after a retention time window is performed in the following steps.
#
##########################################################################################################################################################

	
        resultsPath <- MAIT.object@PhenoData@resultsPath

    
	xsa <- xsAnnotate(rawData(MAIT.object)$xcmsSet)
	xsaF <- groupFWHM(xsa,perfwhm=perfwhm,sigma=sigma)
	cat("Spectrum build after retention time done",fill=TRUE)

##########################################################################################################################################################
#
#  Correlation in the spectra and across samples is calculated and compared to the input correlation threshold value.
#        
##########################################################################################################################################################
        xsaFA <- findIsotopes(xsaF)
	cat("Isotope annotation done",fill=TRUE)

	xsaFA <- groupCorr(xsaFA,cor_eic_th=corrWithSamp,cor_exp_th=corrBetSamp,calcIso=calcIso,calcCiS=calcCiS,calcCaS=calcCaS,pval=pval,graphMethod=graphMethod)
	cat("Spectrum number increased after correlation done",fill=TRUE)

##########################################################################################################################################################
#
#  Once the spectra are finally build, adducts/fragments are obtained (Annotation step).
#        
##########################################################################################################################################################

        if(annotateAdducts==TRUE){
	xsaFA <- findAdducts(xsaFA,rules=adducts,polarity="positive")
	cat("Adduct/fragment annotation done",fill=TRUE)
      }
        
	peakList <- getPeaklist(xsaFA)
	peakList <- peakList[order(as.numeric(peakList$pcgroup)),]
	peakList[,4] <- peakList[,4]/60
	peakList[,5] <- peakList[,5]/60
	peakList[,6] <- peakList[,6]/60
    
	

        peakList[,1:6] <- round(peakList[,1:6],4)
    
##########################################################################################################################################################
#
#  Prints a smaller version of the getPeaklist table where only the mass, the retention time (in minutes) and the spectra index are written.
#        
##########################################################################################################################################################

    
        if(printSpectraTable==TRUE){

          tab<-cbind(peakList[,1],peakList[,4],peakList[,dim(peakList)[2]])
          rownames(tab)<-rep("",dim(tab)[1])
	  colnames(tab)<-c("mass","rt","spectra Index")
			
			if(!file.exists(paste(resultsPath,"Tables",sep="/"))){
		if(resultsPath==""){
                         dir.create("Tables")
                       }else{
			   dir.create(paste(resultsPath,"Tables",sep="/"))
                         }
				
			}else{
				
			   cat(" " ,fill=TRUE)	
			   warning(paste("Warning: Folder",paste(resultsPath,"Tables",sep="/"),"already exists. Possible file overwritting."))
				
			}
	if(resultsPath==""){
	  writeExcelTable(file=tab,file.name="Tables/Spectra")
        }else{
           writeExcelTable(file=tab,file.name=paste(resultsPath,"Tables/Spectra",sep="/"))
         }

          
        }

##########################################################################################################################################################
#
#  The annotation results are saved in the MAIT object which is the output of the function.
#        
##########################################################################################################################################################

    

    xsaFA <- list(xsaFA)
    names(xsaFA) <- "xsaFA"
    
    MAIT.object@RawData@data <- xsaFA

	return(MAIT.object)
  }


