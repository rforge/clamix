\name{empty.symObject}
\alias{empty.symObject}
\title{
Create an Empty Object symObject of Desired Length
}
\description{
Create an empty object \code{symObject} of desired length. Creates vector of lists that represent variables of specified number of categories.
}
\usage{
empty.symObject(nCats)
}
\arguments{
  \item{nCats}{
     vector of number of categories for each variable (length(vector) represents number of variables)
}
}
\details{
}
\value{
Object of type \code{symObject}. Vector of lists that represent variables of specified number of categories.
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
\code{\link{print.symObject}}, \code{\link{summary.symData}}, \code{\link{create.symData}}
}
\examples{
# population tables
object <- empty.symObject(c(18,18))
print(object)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ object }
