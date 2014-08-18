\name{summary.symData}
\alias{summary.symData}
\alias{symdata}
\title{
Summary method for object symData
}
\description{
Summary method for object \code{symData}.
}
\usage{ 
\method{summary}{symData}(object, \dots)  
}
\arguments{
  \item{object}{
     object \code{symData}
}
  \item{\dots}{
}
}
\value{
Object \code{symData}. Represents a list of
  \item{SOs}{vector of symbolic objects (class \code{symObject}). The transformed dataset.}
  \item{so}{An empty \code{symObject}.}
  \item{namedSO}{An empty \code{symObject} with names for categories of each variable.}
  \item{alpha}{A vector of weights}
  \item{type}{A type of \code{symData}. See the argument above.}
}
\author{
Natasa Kejzar
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{create.symData}}, \code{\link{symObject}}
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
\keyword{ methods }
