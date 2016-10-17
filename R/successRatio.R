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


successRatio <- function(classes,tt,ClassWeights){
  
	TP <- 0
	FN <- 0
        SR_class <- vector(length=length(classes))
        
	for (j in c(1:length(classes))){
          
           TP_class<- 0
           FN_class<- 0
        
		for (i in c(1:length(classes))){
                  
			if(i==j){
                          
                                #TP <- (TP+tt[i,i])
				TP <- (TP+ClassWeights[i]*tt[i,i])

                                TP_class <- ClassWeights[i]*tt[i,i]
                                #TP_class <- tt[i,i]
                                
			}else{
                          
				FN <- FN+ClassWeights[i]*tt[i,j]                         
                               # FN_class <- FN
                                FN_class <- FN_class+ClassWeights[i]*tt[i,j]
                                
			 }
                        
                        
		}
                
                SR_class[j] <- TP_class/(FN_class+TP_class)
                
                              
	}
        
	SR <- TP/(FN+TP)

        out <- list(SR,SR_class)
        names(out) <- c("SR","SR_class")
        
	return(out)
        
}
