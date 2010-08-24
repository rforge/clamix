\name{makeEnc}
\alias{makeEnc}
\title{
Function for "Equiprobabilistic" Encoding of Numerical Variables
}
\description{
Prints out "equiprobabilistic" encodings (categorization) for numerical variables.
}
\usage{
makeEnc(var, name, k, file="")
}
\arguments{
  \item{var}{variable values
}
  \item{name}{name of variable
}
  \item{k}{number of categories
}
  \item{file}{file to write encoding into, default "" - to std. output}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The encoding is printed out and is to be used with \code{\link{encodeSO}} function for encoding numerical variables into variables suitable for \code{symObject}.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Vladimir Batagelj
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{encodeSO}}, \code{\link{create.symData}}, \code{\link{symObject}}
}
\examples{
set.seed(42)
makeEnc(runif(100),"Rand",10)
# result - a LIST ():
#encRand <- list(
#  "[0]" = function(x) x<=0,
#  "(0,0.08998052]" = function(x) x<=0.08998052,
#  "(0.08998052,0.2163854]" = function(x) x<=0.2163855,
#  "(0.2163854,0.3590283]" = function(x) x<=0.3590284,
#  "(0.3590283,0.4749971]" = function(x) x<=0.4749972,
#  "(0.4749971,0.6117786]" = function(x) x<=0.6117787,
#  "(0.6117786,0.6932048]" = function(x) x<=0.6932049,
#  "(0.6932048,0.7846928]" = function(x) x<=0.7846929,
#  "(0.7846928,0.914806]" = function(x) x<=0.9148061,
#  "(0.914806,0.9888917]" = function(x) x<=0.9888918,
#  "NA" = function(x) TRUE )
}
\keyword{ manip }
