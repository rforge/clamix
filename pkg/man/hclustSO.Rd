\name{hclustSO}
\alias{hclustSO}
\title{
Adapted hierarchical clustering of modal multi valued modal symbolic data
}
\description{
Function to compute adapted hierarchical clustering for modal multi valued modal symbolic data. Data must be saved in a \code{symData} object.
}
\usage{
hclustSO(dataset,type = "d1")
}
\arguments{
  \item{dataset}{
		\code{symData} object
}
  \item{type}{
  		switch describing dissimilarity measure used "d1" to "d6" or "d1w" (calculating with frequency distributions "fDist") (for measures see Table 1 in reference article below)
}
}
\value{
An object of class \code{\link{hclust}} which describes the tree produced by the clustering process. The object is a list (see hclust) and additional component: 

	\item{leaders}{A vector of \code{symObject}s that represent leaders obtained during merging of clusters}
}
\references{
V. Batagelj, N. Kejzar, and S. Korenjak-Cerne. Clustering of Modal
Valued Symbolic Data. ArXiv e-prints, 1507.06683, July 2015.
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
dataset <- create.symData(datalist,"pDist")
res2006 <- hclustSO(dataset,type="d1") # "d1w" can be used for dataset with "fDist"
plot(res2006,cex=0.2)
rect.hclust(res2006,k=4)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ multivariate }% __ONLY ONE__ keyword per line
\keyword{ cluster }