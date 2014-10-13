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
##     project function is used to get the scores of an spectralData object on the model build by another spectralData object. The needed variables are:
##
##     modelData:   Is a MAIT object whose model's loadings are going to be used to get the scores
##     projectData: Is a MAIT object whose data is going to be projected on the modelData's loadings to get the scores.
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

########################################################################################################################################################## 
##########################################################################################################################################################


project <- function(modelData,
                    projectData){
                     



##########################################################################################################################################################
#
#     The needed data coming from the two spectralData objects is renamed.
#
##########################################################################################################################################################

                    
  
                    data<-scores(projectData)
                    index <- featureID(modelData)
                    idGroup <- featureID(projectData)
                    colnames(data) <- NULL
                    rownames(data) <- NULL

                    
##########################################################################################################################################################
#
#   If index and idGroup are the same, then index is redifined to be the positions where the idGroup changes its value
#
##########################################################################################################################################################

                    
                       idGroup<-as.numeric(idGroup)
                    
                       if(length(index)==length(idGroup)&&sum(index-idGroup)==0){
                        
                          index<-idGroup[which(diff(idGroup)!=0)]
                         
                          if(idGroup[length(idGroup)]!=idGroup[length(idGroup)-1]){
                           
                             index<-c(index,idGroup[length(idGroup)])
                            
                          }
                       }

                    
##########################################################################################################################################################
#
#   If no dimensionality reduction method is used, the scores are just the variables contained in index
#
##########################################################################################################################################################

                    
                       
                       out <- data[index,]
                       
                    

         
       
                    return(out)
                    
                  }
