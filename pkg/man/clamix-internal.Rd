\name{clamix-internal}
\alias{acceptbin}
\alias{acceptinterval}
\alias{distCl}
\alias{distSO}
\alias{testSOvar}
\alias{testSOvarP}
\alias{testSOvarQ}
\title{
Internal Functions in the clamix package}
\description{Internal functions of the package.}
\usage{
acceptbin(x,n,p)
acceptinterval(x, n, level, tolerance=1e-04)
distSO(X,Y,nVar,alpha)
distCl(X,Y,nVar,alpha)
testSOvar(SO, var, total, a)
testSOvarP(SO, var, total, a)
testSOvarQ(SO, var, total, a, qmin)
}
\arguments{
  \item{x}{observed unit of Binomial distribution}
  \item{n}{number of trials in Bin(n,p)}
  \item{p}{probability in Bin(n,p)}
  \item{level}{= (1-alpha) interval for p in Bin(n,p)}
  \item{tolerance}{tolerance for computing the interval}
  \item{X}{unit - \code{symObject}}
  \item{Y}{unit - \code{symObject}}
  \item{nVar}{number of variables}
  \item{alpha}{weight vector}
  \item{SO}{unit (leader) of class \code{symObject}}
  \item{var}{number of variable to check}
  \item{total}{unit made of summed total dataset}
  \item{a}{alpha value for acceptance interval}
  \item{qmin}{the minimal size of quotient between acceptance probability and actual probability}
}
\keyword{ internal }
