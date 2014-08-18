\name{create.symData}
\alias{create.symData}
\title{
Creates object symData from a list of data.frames
}
\description{
Creates object \code{symData} from a named list of data.frames. Each should \code{data.frame} represent one variable of the units.
}
\usage{
create.symData(datalist, type = "gDist", alpha = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datalist}{
     a named list of \code{data.frame}s. Each \code{data.frame} should represent one variable. Data.frames should be named as variables. Each dataframe should have the same number of rows (same number of units).
}
  \item{type}{
     type of \code{symData}. \dQuote{pDist} for probability, \dQuote{fDist} for frequency and \dQuote{gDist} for general data distribution.
}
  \item{alpha}{
     a vector of weights for each variable (length(alpha) should be the same as number of data.frames). Default \code{NULL} - variables have equal weights.
}
}
\value{
Object \code{symData}. Represents a list of
  \item{SOs}{vector of symbolic objects (class \code{symObject}). The transformed dataset.}
  \item{so}{An empty \code{symObject}.}
  \item{namedSO}{An empty \code{symObject} with names for categories of each variable.}
  \item{alpha}{A vector of weights}
  \item{type}{A type of \code{symData}. See the argument above.}
%% ...
}
\author{
   Natasa Kejzar
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{summary.symData}}, \code{\link{print.symObject}}, \code{\link{empty.symObject}}
}
\examples{
data(popul06f)
data(popul06m)
datalist <- list("M"=popul06f,"F"=popul06m)
dataset <- create.symData(datalist,"fDist")
summary(dataset)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }