\name{hclustSO}
\alias{hclustSO}
\title{
Adapted hierarchical clustering of modal multi valued modal symbolic data
}
\description{
Function to compute adapted hierarchical clustering for modal multi valued modal symbolic data. Data must be saved in a \code{symData} object.
}
\usage{
hclustSO(dataset,type = "euclid")
}
\arguments{
  \item{dataset}{
		\code{symData} object
}
  \item{type}{
  		switch describing dissimilarity measure used (currently "euclid" and "d4" (see Table 1 in reference article below) implemented)
}
}
\value{
An object of class \code{\link{hclust}} which describes the tree produced by the clustering process. The object is a list (see hclust) and additional component: 

	\item{leaders}{A vector of \code{symObject}s that represent leaders obtained during merging of clusters}
}
\references{
Kejzar,~N., Korenjak-Cerne,~S., and Batagelj,~V., Clustering of distributions: A case of patent citations, Journal of Classification, 28, 2, p. 156--183, (2011)
}
\author{
	Vladimir Batagelj
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