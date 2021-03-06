\name{spectralSigFeatures}
\alias{spectralSigFeatures}
\title{Extract significant features from a MAIT object}
\description{
Function spectralSigFeatures takes a \code{\link{MAIT-class}} object and obtains which of the variables are significant given a p-value threshold. The parameters of the significant features can ve printed to an output table (TRUE by default). Depending on the number of classes in the data, the function chooses between using ANOVA tests through function \link{spectralAnova}, or T-Student tests by using function \link{spectralTStudent}.
}
\usage{

spectralSigFeatures(MAIT.object = NULL,
                    pvalue = 0.05,
                    p.adj = "none", 
                    printCSVfile = FALSE,
                    scale = FALSE,
                    parametric = TRUE,
                    var.equal = FALSE, 
                    test.fun = NULL,
                    jitter = FALSE,
                    jitter.factor = 1,
                    jitter.amount = 0, 
                    namefun = NULL)
}
\arguments{
  \item{MAIT.object}{
A \link{MAIT-class} object where function \link{peakAggregation} has already been applied. The output of the function is going to be an update of the same \link{MAIT-class} object.
}
  \item{pvalue}{
P-value threshold. Variables having a p-value lower than this value is considered as a significant variable.
}
  \item{p.adj}{
Character with the name of the posthoc method to be applied to correct the pvalues. The supported methods are that of the \link{p.adjust} function
}

  \item{printCSVfile}{
Set to TRUE if an output table has to be produced. See function \link{sigPeaksTable} for more information.
}
  \item{scale}{
Set to FALSE by default. When set to TRUE, a unit variance scaling of the data when no peak aggregation is performed. If a peak aggregation method is applied, this parameter is ignored.
}

  \item{parametric}{
If it is set to TRUE, the statistical tests to be applied will be parametrical tests (e.g. ANOVA, TStudent or Welch's tests). Non-parametrical tests (e.g. Kruskal-Wallis, Mann-Whitney tests) are applied otherwise.
}

  \item{var.equal}{
Set to FALSE by default. When set to TRUE, a Student's T-Test is applyied when having 2 classes in the data. If it is set to FALSE, a Welch's test is applyied instead.
}

  \item{test.fun}{
Function of the user-defined posthoc method to be applyied.
}

  \item{jitter}{
If it is set to TRUE, a jitter noise is added to the data. This is useful when applying Mann-Whitney tests with ties.
}

  \item{jitter.factor}{
See argument factor of the function \link{jitter}.
}

  \item{jitter.amount}{
See argument amount of the function \link{jitter}.
}


  \item{namefun}{
Name of the user-defined posthoc test in the argument test.fun.
}


}
\value{
  A \link{MAIT-class} object containing the significant features of the scores slot of \link{MAIT-class} object used as an input.
}
\seealso{
  \code{\link{MAIT-class}} \code{\link{spectralTStudent}} \code{\link{spectralAnova}} \code{\link{sigPeaksTable}}
}

\examples{
data(MAIT_sample)
MAIT<-spectralSigFeatures(MAIT,p.adj="fdr",parametric=TRUE)
}

\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}
\keyword{file}
