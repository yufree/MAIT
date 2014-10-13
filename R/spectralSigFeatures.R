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


spectralSigFeatures <- function(MAIT.object=NULL,
                           pvalue=0.05,
                           p.adj="none",
                           printCSVfile=FALSE,
                           scale=FALSE,
                           parametric=TRUE,
                           var.equal=FALSE,
                           test.fun=NULL,
                           jitter=FALSE,
                           jitter.factor=1,
                           jitter.amount=0,
                           namefun=NULL)
{

  if(length(classes(MAIT.object))==1){
    if(is.na(classes(MAIT.object))==TRUE){

  stop("No class information saved in the MAIT object.")
}
}
  
   if(method(MAIT.object)==""){
      cat("Skipping peak aggregation step...")
      MAIT.object <- peakAggregation(MAIT.object,scale=scale)
   }

if(is.null(test.fun)){

   
   if(length(classes(MAIT.object))==2){
      if(parametric==TRUE){
         if(var.equal==TRUE){
           
            out<-spectralTStudent(pvalue=pvalue,
                           p.adj=p.adj,
                           MAIT.object=MAIT.object,
                           printCSVfile=printCSVfile)
            
         }else{
              
            out<-spectralWelch(pvalue=pvalue,
                           p.adj=p.adj,
                           MAIT.object=MAIT.object,
                           printCSVfile=printCSVfile)
         }
         
      }else{
        
         out<-spectralWilcox(pvalue=pvalue,
                           p.adj=p.adj,
                           MAIT.object=MAIT.object,
                           printCSVfile=printCSVfile,
                           jitter=jitter,
                           jitter.factor=jitter.factor,
                           jitter.amount=jitter.amount)
      }

    
   }else{

      if(parametric==TRUE){
        
         out<-spectralAnova(pvalue=pvalue,
                           p.adj=p.adj,
                           MAIT.object=MAIT.object,
                           printCSVfile=printCSVfile)
         
      }else{
        
      out<-spectralKruskal(pvalue=pvalue,
                           p.adj=p.adj,
                           MAIT.object=MAIT.object,
                           printCSVfile=printCSVfile)
      }
   }
 

}else{

   out<-spectralFUN(pvalue=pvalue,
                           p.adj=p.adj,
                           MAIT.object=MAIT.object,
                           printCSVfile=printCSVfile,
                           test.fun=test.fun,
                           namefun=namefun)
   }

   if(length(featureSigID(out))==0){

      warning("No significative features found with the selected parameters.")
  
   }else{

      aux <- sigPeaksTable(out,printCSVfile=printCSVfile)
      
   }
   
   
return(out)

}
