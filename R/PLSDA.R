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


PLSDA <- function (Xtrain, Ytrain, Xtest = NULL, ncomp, nruncv = 0, alpha = 2/3,
priors = NULL) {
    ntrain <- nrow(Xtrain)
    Ytrain <- as.factor(Ytrain)
    if (is.vector(Xtest)) {
        Xtest <- matrix(Xtest, 1, length(Xtest))
    }
    if (is.null(Xtest)) {
        Xtest <- Xtrain
    }
    if (nruncv == 0 & length(ncomp) > 1)
	stop("Since length(ncomp)>1, nruncv must be >0")
    if (nruncv > 0) {
        ncomp <- pls.lda.cv(Xtrain, Ytrain, ncomp = ncomp, nruncv = nruncv,
							alpha = alpha, priors = priors)
    }
    pls.out <- pls.regression(Xtrain = Xtrain, Ytrain = transformy(Ytrain),
							  Xtest = NULL, ncomp = ncomp)
    Ztrain <- as.data.frame(matrix(pls.out[[4]], ntrain, ncomp))
    Ztrain$y <- Ytrain
    Ztest <- as.data.frame(scale(Xtest, center = pls.out$meanX,
								 scale = FALSE) %*% pls.out$R)
    if (is.null(priors)) {
        lda.out <- lda(formula = y ~ ., data = Ztrain)
    }
    else {
        lda.out <- lda(formula = y ~ ., data = Ztrain, prior = priors)
    }
    predclass <- predict(object = lda.out, newdata = Ztest)
    return(list(predclass = predclass, ncomp = ncomp,regressionPLS=pls.out))
}
