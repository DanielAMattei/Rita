% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/master_Functions.R
\name{MasterXform}
\alias{MasterXform}
\title{Master Transformation Function}
\usage{
MasterXform(c, data)
}
\arguments{
\item{c}{Input specifying the test to run (scalar)}

\item{data}{The data of a univariate distribution for which the test statistic is computed (vector)}
}
\value{
Output from the appropriate subfunction (list)
}
\description{
This is a master function used to perform the appropriate transformation(s) within
the 'Rita' function.
}
\examples{
values <- rnorm(100)
x <- MasterXform(c = 2, data = values)
}
