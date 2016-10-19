\name{clamix-package}
\alias{clamix-package}
\alias{clamix}
\docType{package}
\title{
Clustering modal valued symbolic data
}
\description{
Two clustering methods for symbolic objects are implemented: the adapted leaders method and the adapted agglomerative hierarchical clustering Ward's method. Additional: classes \code{symObject} and \code{symData} for saving symbolic objects with \code{print} and \code{summary} (for symData) methods and functions to help producing such objects are implemented.

Additional functions to help cluster descriptions that are implemented: \code{computeTotal}, \code{create.symData}, \code{empty.symObject}, \code{encodeSO}, \code{makeEnc}, \code{specContrasts}.

Datasets: \code{popul06f}, \code{popul06m}, a list of data.frames: \code{TIMSSSlo}. 
}
\details{
\tabular{ll}{
Package: \tab clamix\cr
Type: \tab Package\cr
Version: \tab 1.02\cr
Date: \tab 2016-10-19\cr
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
V. Batagelj, N. Kejzar, and S. Korenjak-Cerne. Clustering of Modal
Valued Symbolic Data. ArXiv e-prints, 1507.06683, July 2015.}
\keyword{ package }
%%\seealso{
%%\code{\link[clustDDist:clustDDist-package]{clustDDist}}
%%}
