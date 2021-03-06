\name{encodeSO}
\alias{encodeSO}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Encode numerical vector using given encoding
}
\description{
Function encodes numerical vector using given encoding. Encoding could be obtained by function \code{makeEnc}.
}
\usage{
encodeSO(x, encoding, codeNA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{number    
}
  \item{encoding}{function name to produce the encoding. See result of \code{makeEnc}.
}
  \item{codeNA}{code (or number) for \code{NA}.
}
}
\value{
Returns category into which the number \code{x} is encoded.}
\author{
Vladimir Batagelj
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{makeEnc}
}
\examples{
set.seed(42)
testset <- runif(100)
makeEnc(testset,"Rand",10,file="temp.R") # make "encRand" encoding
source("temp.R")
unlink("temp.R") # tidy up
testcat <- sapply(testset,function(x) encodeSO(x,encRand,NA))
tabulate(testcat)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
