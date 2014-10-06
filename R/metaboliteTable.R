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


metaboliteTable<-function(MAIT.object,printCSVfile=FALSE){
  if(printCSVfile==TRUE){    
                                    	
	   if(!file.exists(paste(resultsPath(MAIT.object),"Tables",sep="/"))){
		
		   dir.create(paste(resultsPath(MAIT.object),"Tables",sep="/"))
                   writeExcelTable(metaboliteTable(MAIT.object),paste(resultsPath(MAIT.object),"Tables","SearchTable",sep="/"))

		
	   }else{
		
		cat(" " ,fill=TRUE)
		warning(paste("Folder",paste(resultsPath(MAIT.object),"Tables",sep="/"),"already exists. Possible file overwritting.",sep=" "))
          	write.table(as.data.frame(metaboliteTable(MAIT.object)),paste(paste(resultsPath(MAIT.object),"Tables","SearchTable",sep="/"),".csv",sep=""),col.names=TRUE,row.names=FALSE,sep=",")

		
	}
     }
return(MAIT.object@FeatureInfo@metaboliteTable)}
