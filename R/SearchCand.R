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

 
SearchCand <- function(candidate,
                       dataBase,
                       peakTolerance){

	dataBase <- as.matrix(dataBase[order(as.numeric(as.character(dataBase[,4]))),])
	aux <- c(as.numeric(dataBase[,4]),(candidate-peakTolerance))
#	lowIndex<- which(order(as.numeric(aux))==length(aux))
        lowIndex <- which(sort(aux)%in%(candidate-peakTolerance))
	aux <- c(as.numeric(dataBase[,4]),candidate+peakTolerance)
        highIndex <- which(sort(aux)%in%(candidate+peakTolerance))
	#highIndex<- which(order(as.numeric(aux))==length(aux))
	if(lowIndex!=highIndex){
          
		unk <- 0
                
		if(!is.matrix(dataBase[lowIndex:(highIndex-1),])){
                  
			dataBase[lowIndex:(highIndex-1),] <- matrix(dataBase[lowIndex:(highIndex-1),],nrow=1)
                        
		}
                
		out <- list(dataBase[lowIndex:(highIndex-1),],unk)
		names(out) <- c("SearchCand","unk")
                
	}else{
          
		unk <- 1
		out <- list(as.matrix(rep("unknown",5),nrow=1),unk)
		names(out) <- c("SearchCand","unk")
                
	}
        
	return(out)
        
}
