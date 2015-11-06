\name{TIMSSSlo}
\alias{TIMSSSlo}
\docType{data}
\title{
Data from TIMSS (Trends in International Mathematics and Science Study Advanced 2008) questionnaire from 2008 for Mathematics for country Slovenia. The data are about teacher and students Math background.
}
\description{
A \code{list} of \code{data.frame}s. Statistical units are teachers. Each dataframe represents one variable from the questionnaire, "variables" represent this variable categories. There are 147 variables for teachers and 77 variables for students (taught by that teacher). The names of variables are the the names of list items and are thoroughly explained in the TIMSS Advanced 2008 documentation at the website of IEA.

The variables are mostly categorical, the numerical ones were categorized according to the 5 equidistant quantiles. One teacher and 20 students were excluded due to many missing values.
}
\usage{data("TIMSSSlo")}
\format{
  A list of 224 data.frames (variables). The first data.frame is thus described with name "MT2GAGE"
  \describe{
    \item{\samp{UNDER 25}}{a numeric vector}
    \item{\samp{25-29}}{a numeric vector}
    \item{\samp{30-39}}{a numeric vector}
    \item{\samp{40-49}}{a numeric vector}
    \item{\samp{50-59}}{a numeric vector}
    \item{\samp{60 OR OLDER}}{a numeric vector}
    \item{\samp{NA}}{a numeric vector}
  }
}
\source{
The International Association for the Evaluation of Educational Achievement (IEA) \url{http://www.iea.nl/data.html}
}
\examples{
data(TIMSSSlo)
names(TIMSSSlo)
SymbolicData = create.symData(TIMSSSlo)
}
\keyword{datasets}
