\name{as.symData}
\alias{as.symData}
\title{
Creates object symData from a list if possible
}
\description{
Creates object \code{symData} from a list of \code{symObject}s. If not returns \code{FALSE}.
}
\usage{
as.symData(L, nvnames = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{L}{
     a list of \code{symObject}s.
}
  \item{nvnames}{
     names of components of each variable in \code{symObject}.
}
}
\value{
Either \code{FALSE} or object \code{symData}. Object represents a list of
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
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
