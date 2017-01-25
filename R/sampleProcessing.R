########################################################|
##   Metabolite Automatic Identifimessageion Toolkit (MAIT) |
##                                                      |
##                                                      |
##   written by Francesc Fernández Albert               |
##   contact mail: francesc.fernandez.albert@upc.edu    |
##   date: 7/29/2013                                    |
##                                                      |
##   SISBIO Group. ESAII Department                     |
##   Technical University of messagealonia                  |
##   Nutrition Department                               |
##   University of Barcelona                            |
##                                                      |
##   SISBIO Group. ESAII Department                     |
##   Universitat Politècnica de messagealunya               |
##   ___________________________________________________|
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
#'   sampleProcessing function takes a set of netCDF files containing LC/MS sample data and performs a peak detection, retention time correction and peak grouping steps using the xcms package. A MAIT object is created and all the informated is saved in it. Each of these steps are done using the package xcms.
#' @param dataDir The netCDF sample files of each class present in the data should be stored in a folder called /(working directory)/Data/(ClassName) replacing (ClassName) for the name of the folder where the files are stored.
#' @param retcorrMethod Method used to correct the retention times values of the variables. By default is set to "loess".
#' @param groupMethod  Method used to build the group peaks of variables. By default is set to "density".
#' @param filterMethod Filtering method applied in the peak detection step. (Set to "matchedFilter" by default).
#' @param project Project folder name under which the results will be saved. This folder will be created in the working directory.

sampleProcessing <- function(dataDir = NULL,
                             project = NULL,
                             index = F,
                             BPPARAM=SnowParam(workers = 12),
                             xsmethod = "centWave",
                             peakwidth = c(14, 25),
                             ppm = 2.5,
                             noise = 0,
                             snthresh = 10,
                             mzdiff = -0.00395,
                             prefilter = c(3, 100),
                             mzCenterFun = "wMean",
                             integrate = 1,
                             fitgauss = FALSE,
                             verbose.columns = FALSE,
                             rmethod = "obiwarp",
                             plottype = "none",
                             distFunc = "cor_opt",
                             profStep = 1,
                             center = 2,
                             response = 1,
                             gapInit = 0.6176,
                             gapExtend = 2.4,
                             factorDiag = 2,
                             factorGap = 1,
                             localAlignment = 0,
                             gmethod = "density",
                             bw = 0.25,
                             mzwid = 0.0021748,
                             minfrac = 1,
                             minsamp = 1,
                             gmax = 50,
                             ...)
{
        if (is.null(dataDir)) {
                stop("No input directory was given")
        }
        
        if (is.null(project)) {
                stop("No project name was included")
        }
        
        
        
        parameters <- list(dataDir,
                           project,
                           xsmethod,
                           peakwidth,
                           ppm,
                           prefilter,
                           rmethod,
                           bw,
                           mzwid)
        
        names(parameters) <- c(
                "dataDir",
                "project",
                "filterMethod",
                "centWave peakwidth",
                "ppm",
                "prefilter",
                "retcorrMethod",
                "bw",
                "mzwid"
        )
        
        MAIT.object <- new("MAIT")
        
        
        MAIT.object@RawData@parameters@sampleProcessing <-
                parameters
        writeParameterTable(parameters(MAIT.object), folder = resultsPath(MAIT.object))
        
        class <-  list.files(dataDir)
        classNum <- vector(length = length(class))
        fileList <-
                list.files(path = paste(dataDir, list.files(path = dataDir), sep = "/"),
                           full.names = TRUE)
        
        for (i in 1:length(class)) {
                classNum[i] <-
                        length(list.files(paste(
                                dataDir, list.files(dataDir)[i], sep = "/"
                        )))
        }
        
        
        classes <- rep(class, classNum)
        
        if (length(list.files(dataDir)) == 1) {
                warning("Warning: Input data only has one class!")
        }
        
        if (is.null(project)) {
                warning("Warning: Project name is empty!")
        }
        #  The results folder is created if it does not exist
        if (!is.null(project)) {
                resultsPath <- paste("Results", project, sep = "_")
                dir.create(resultsPath)
                
        } else{
                resultsPath <- "Results"
                dir.create(resultsPath)
                
        }
        #  Peak detection step.
        peaks <-
                xcmsSet(
                        files = fileList,
                        snthresh = snthresh,
                        method = xsmethod,
                        ppm = ppm,
                        peakwidth = peakwidth,
                        method = xsmethod,
                        snthresh = snthresh,
                        mzdiff = mzdiff,
                        BPPARAM = BPPARAM,
                        noise = noise,
                        prefilter = prefilter,
                        mzCenterFun = mzCenterFun,
                        integrate = integrate,
                        fitgauss = fitgauss,
                        verbose.columns = verbose.columns,...
                )
        
        message("Peak detection done", fill = TRUE)
        #  Peak grouping and retention time correction steps.
        
        groups <- group(
                peaks,
                method = gmethod,
                bw = bw,
                mzwid = mzwid,
                minfrac = minfrac,
                minsamp = minsamp,
                max = gmax
        )
        
        if (rmethod != "none") {
                retcorr_groups <-
                        retcor(
                                groups,
                                method = rmethod,
                                plottype = plottype,
                                distFunc = distFunc,
                                profStep = profStep,
                                center = center,
                                response = response,
                                gapInit = gapInit,
                                gapExtend = gapExtend,
                                factorDiag = factorDiag,
                                factorGap = factorGap,
                                localAlignment = localAlignment
                        )
                message("Retention time correction done", fill = TRUE)
                groups <-
                        group(
                                retcorr_groups,
                                method = gmethod,
                                bw = bw,
                                mzwid = mzwid,
                                minfrac = minfrac,
                                minsamp = minsamp,
                                max = gmax
                        )
                message("Peak grouping after samples done", fill = TRUE)
                
        } else{
                message("Skipping retention time correction...", fill = TRUE)
                
        }
        #  Fill missing peaks step
        
        fPeaks <- fillPeaks(groups,BPPARAM = BPPARAM)
        message("Missing Peak integration done", fill = TRUE)
        
        #  Writting the results in a new MAIT object
        fPeaks <- list(fPeaks)
        names(fPeaks) <- "xcmsSet"
        
        MAIT.object@RawData@data <- fPeaks
        MAIT.object@PhenoData@classes <- class
        MAIT.object@PhenoData@classNum <- classNum
        MAIT.object@PhenoData@resultsPath <- resultsPath
        
        return(MAIT.object)
}
