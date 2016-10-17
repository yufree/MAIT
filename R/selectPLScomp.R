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


selectPLScomp <- function(data, class, max.comp) {
	error.comp <- rep(-9, max.comp)
	for (i in 1:max.comp) {
		error.comp[i] <- length(which(class != PLSDA(Xtrain=data, Ytrain=class, ncomp = i)$predclass$class))
	}
	return(which.min(error.comp))
}
