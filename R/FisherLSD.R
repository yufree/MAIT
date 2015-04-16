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


FisherLSD <- function(data,classes,index,DFerror,MSerror,numClasses){
	LSD <- LSD.test(y=data[index,],trt=classes,DFerror=DFerror,MSerror=MSerror);
	LSD <- LSD$groups[order(LSD$groups$trt),]
	groups <- paste(LSD$M[1],LSD$M[2],sep=" ")
	for (i in c(3:(numClasses))){
		groups <- paste(groups,LSD$M[i],sep=" ")
	}
	out <- list(LSD$trt,groups,LSD$means)
	names(out) <- c("trt","group","means")
	return(out)
}
