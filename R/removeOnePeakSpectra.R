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


removeOnePeakSpectra <- function(data,
                                 idGroup){


########################################################################################################################################################## 
##########################################################################################################################################################
##       
##   removeOnePeakSpectra function takes a numerical matrix along with a vector having the variables' ids and removes those variables whose id is unique (i.e. not shared   ##   with more variables . The input variables are:
##   
##
##   data:          Numerical matrix or dataset containing the peak intensity.
##   idGroup:       Numeric vector containing the spectral ID of each peak present in data.
##   signVariables: If it is a numeric vector, only the variables with the ID number present in the vector are chosen (set to NULL by default).
##
########################################################################################################################################################## 
##########################################################################################################################################################



  
        data <- as.data.frame(data)
        index<-idGroup[which(diff(idGroup)!=0)]
        spectra <- matrix(nrow=1,ncol=ncol(data))
        peaks <- vector(length=1)
        
        if(idGroup[length(idGroup)]!=idGroup[length(idGroup)-1]){
          
           index<-c(index,idGroup[length(idGroup)])
           
        }
  
        for (i in c(1:length(index))){


##########################################################################################################################################################
#
#   Variable spectrum contains the peak matrix of the spectrum whose id number is equal to index[i]
#                
##########################################################################################################################################################

          
           spec <- as.matrix(data[which(index[i]==idGroup),])



##########################################################################################################################################################
#
#   Only those spectra having more than one peak are attached to the output matrix.
#                
##########################################################################################################################################################


          
           if(dim(spec)[1]>1){
             
              spectra <- rbind(spectra,spec)
              peaks <- c(peaks,rep(index[i],length(which(index[i]==idGroup))))
              
           }
           
        }
        
        spectra <- spectra[-1,]
        peaks <- peaks[-1]
        out <- list(spectra,peaks)
        names(out) <- c("spectra","idGroup")
        return(out)
        
}
