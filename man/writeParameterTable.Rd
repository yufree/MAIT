\name{writeParameterTable}
\alias{writeParameterTable}
\title{Writes a csv table containing the parameters launched in the MAIT analysis}
\description{
Function \code{writeParameterTable} writes a csv table where all the provided input parameters in the whole MAIT analysis are saved.
}
\usage{

writeParameterTable(listParameters,
         folder)

}

\arguments{
  \item{listParameters}{
The list of parameters to be printed. This input should be an object of the class MAIT.Parameters
}

  \item{folder}{
The folder where the csv file is going to be saved
}
}


\value{
A csv file containg the input parameters of the whole run.
}

\author{Francesc Fernandez,
 \email{francesc.fernandez.albert@upc.edu}}
\keyword{internal}