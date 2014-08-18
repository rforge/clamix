\name{leaderSO}
\alias{leaderSO}
\title{
Adapted leaders clustering of modal multi valued symbolic data
}
\description{
Function to compute adapted leaders clustering for modal multi valued symbolic data. Data must be saved in a \code{symData} object.
}
\usage{
leaderSO(dataset, maxL, initial = NULL, stabil = 1e-06, report = TRUE, interact = TRUE,type = "euclid")
}
\arguments{
  \item{dataset}{
		\code{symData} object
}
  \item{maxL}{
		number of clusters to produce
}
  \item{initial}{
		provide initial leaders for clustering; if NULL they are selected randomly
}
  \item{stabil}{
		stability of the results (how close should be the minimization functions of two algorithm steps to exit)
}
  \item{report}{
		produce report while working
}
  \item{interact}{
		do you want interaction during algorithm steps
}
  \item{type}{
  		switch describing dissimilarity measure used (currently "euclid" and "d4" (see Table 1 in reference article below) implemented)
}
}
\value{
Returned is a list of elements:
list(clust=clust,leaders=L,R=R,p=p)
	\item{clust}{partition of clusterings}
	\item{leaders_symData}{\code{symData} object of leaders for each cluster (for possible direct input into \code{hclustSO} function)}
	\item{R}{minimal error (distance to leader) for each cluster}
	\item{p}{sum of errors for each cluster}
}
\references{
Kejzar,~N., Korenjak-Cerne,~S., and Batagelj,~V., Clustering of distributions: A case of patent citations, Journal of Classification, 28, 2, p. 156--183, (2011)
}
\author{
	Vladimir Batagelj
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{hclustSO}}, \code{\link{create.symData}},\code{\link{print.symObject}}
}
\examples{
data(popul06f)
data(popul06m)
datalist <- list("M"=popul06f,"F"=popul06m)
dataset <- create.symData(datalist,"fDist")
\dontrun{
res2006 <- leaderSO(dataset,5)
# Times repeat (if <=0, no repetitions)
summary(res2006$leaders_symData)
print(res2006$leaders_symData$SOs)
}
}
\keyword{ models }
\keyword{ multivarate }% __ONLY ONE__ keyword per line
\keyword{ cluster }