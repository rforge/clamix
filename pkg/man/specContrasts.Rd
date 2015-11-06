\name{specContrasts}
\alias{specContrasts}
\title{
Compute specificities (importances) and contrasts ("partial" importances) between two clusters.
}
\description{
Computes specificities for a variable comparing two cluster leaders and contrasts (partial importances) of the second cluster leader towards first for each variable component.
}
\usage{
specContrasts(rs,rc,j)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rs}{
     first cluster leader in the form of \code{symObject}, usually leader of whole dataset.
}
  \item{rc}{
     second cluster leader in the form of \code{symObject}.
}
  \item{j}{
     an index (or a name) of a variable.
}
}
\value{
A list of
  \item{specificity}{A value.}
  \item{contrasts}{A vector of components of NaN, Inf (see reference) or a value x where |x| > 1.}
%% ...
}
\references{
Batagelj,~V., Kejzar,~N. and Korenjak-Cerne,~S., Clustering of Modal Valued Data, a chapter in Analysis of Distributional Data (Ed. Brito,~P.), (submitted)
}
\author{
   Natasa Kejzar
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{hclust}}, \code{\link{hclustSO}}, \code{\link{leaderSO}}, \code{\link{print.symObject}}, \code{\link{empty.symObject}}
}
\examples{
data(TIMSSSlo)
tsdata = create.symData(TIMSSSlo)
dendro <- hclustSO(tsdata)
rs = dendro$leaders[[nrow(dendro$merge)]] # the whole dataset
rc = dendro$leaders[[84]] # first cluster from reference chapter

specContrasts(rs,rc,1) # first variable ("MT2GTAUT") 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }