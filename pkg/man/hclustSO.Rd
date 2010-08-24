\name{hclustSO}
\alias{hclustSO}
\title{
Adapted Hierarchical Clustering of Symbolic Data Described as Discrete Distributions
}
\description{
Function to compute adapted hierarchical clustering for symbolic data described by discrete distributions. Data must be saved in a \code{symData} object.
}
\usage{
hclustSO(dataset)
}
\arguments{
  \item{dataset}{
		\code{symData} object
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
An object of class \code{\link{hclust}} which describes the tree produced by the clustering process. The object is a list (see hclust) and additional component: 

	\item{leaders}{A vector of \code{symObject}s that represent leaders obtained during merging of clusters}
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
\code{\link{leaderSO}}, \code{\link{create.symData}}, \code{\link{hclust}}, \code{\link{print.symObject}}
}
\examples{
data(popul06f)
data(popul06m)
datalist <- list("M"=popul06f,"F"=popul06m)
dataset <- create.symData(datalist,"fDist")
res2006 <- hclustSO(dataset)
plot(res2006,cex=0.2)
rect.hclust(res2006,k=4)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ multivariate }% __ONLY ONE__ keyword per line
\keyword{ cluster }