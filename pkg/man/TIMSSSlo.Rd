\name{TIMSSSlo}
\alias{TIMSSSlo}
\docType{data}
\title{
Data from TIMSS (Trends in International Mathematics and Science Study Advanced 2008) questionnaire from 2008 for Mathematics for country Slovenia. The data are about teacher (93) and students (2120) Math background.
}
\description{
A \code{list} of \code{data.frame} s. Statistical units are teachers. Each dataframe represents one variable from the questionnaire, "variables" represent this variable categories. There are 2 variables for teachers and 3 variables for students (taught by that teacher). The names of variables are the the names of list items and are thoroughly explained in the TIMSS Advanced 2008 documentation at the website of IEA.

The variables are mostly categorical, the numerical ones were categorized according to the 5 equidistant quantiles. One teacher and 20 students were excluded due to many missing values.
}
\usage{data("TIMSSSlo")}
\format{
  A list of 5 data.frames (variables). The variables (data.frames) are 
  \describe{
  \item{"MT2TAUT"}{teachers' years teaching: [2,10], (10,15], (15,21], (21,26.6],(26.6,35]}
  \item{"MT2GOTDC"}{teachers' professional interaction with coworkers: never or almost never, 2-3 times per month, 1-3 times per week, daily or almost daily}
  \item{"MS2GSEX"}{students' gender: GIRL, BOY}
  \item{"MS2MACRH"}{quantification of frequency of reviewing homework: every or almost every lesson, about half the lessons, some lessons, never}
  \item{"MS2MHTIM"}{minutes per week of doing homework: [0,30], (30,60], (60,120], (120,180], (180,600]}
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
