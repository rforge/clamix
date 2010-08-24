\name{print.symObject}
\alias{print.symObject}
\alias{symObject}
\title{
Print Method for Object symObject
}
\description{
Print method for object \code{symObject}
}
\usage{ 
\method{print}{symObject}(x, \dots)  
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     object \code{symObject}
}
  \item{\dots}{
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Object \code{symData}. Vector of lists that represent variables of specified number of categories.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Natasa Kejzar
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{empty.symObject}}, ~~~
}
\examples{
object <- empty.symObject(c(18,18))
print(object)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
