\name{clamix-package}
\alias{clamix-package}
\alias{clamix}
\docType{package}
\title{
Clustering modal valued symbolic data
}
\description{
Two clustering methods for symbolic objects are implemented: the adapted leaders method and the adapted agglomerative hierarchical clustering Ward's method. Additional: classes \code{symObject} and \code{symData} for saving symbolic objects with \code{print} and \code{summary} (for symData) methods and functions to help producing such objects are implemented.
}
\details{
\tabular{ll}{
Package: \tab clamix\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-08-14\cr
License: \tab GPL-2\cr
}
Adapted agglomerative and hierarchical clustering methods are compatible - they can be viewed as two approaches for solving the same clustering optimization problem. In the obtained clustering to each cluster is assigned its leader. The descriptions of the leaders offer simple interpretation of the clusters' characteristics. The leaders method enables us to efficiently solve clustering problems with large number of units; while the agglomerative method is applied on the obtained leaders and enables us to decide upon the right number of clusters on the basis of the corresponding dendrogram.
}
\author{
Vladimir Batagelj

Maintainer: Natasa Kejzar <natasa.kejzar@mf.uni.lj.si>
Vladimir Batagelj <vladimir.batagelj@fmf.uni-lj.si>
}
\references{
~~ Literature or other references for background information ~~
}
\keyword{ package }
%%\seealso{
%%\code{\link[clustDDist:clustDDist-package]{clustDDist}}
%%}
