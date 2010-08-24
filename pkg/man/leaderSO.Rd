\name{leaderSO}
\alias{leaderSO}
\title{
Adapted Leaders Clustering of Symbolic Data Described as Discrete Distributions
}
\description{
Function to compute adapted leaders clustering for symbolic data described by discrete distributions. Data must be saved in a \code{symData} object.
}
\usage{
leaderSO(dataset, maxL)
}
\arguments{
  \item{dataset}{
		\code{symData} object
}
  \item{maxL}{
		number of clusters to produce
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
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