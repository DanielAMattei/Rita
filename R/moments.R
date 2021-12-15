
#MOMENTS FUNCTIONS##################################################################################
#' Adjusted Fisher-Pearson Skewness Coefficient with Sample-size Correction Factor
#'
#' @param data The data for which skewness is computed (vector)
#' @param sd The population standard deviation, used to compute skewness (scalar)
#'
#' @return The skewness value (scalar)
#' @export
#'
#' @references Shreve, Joni N. and Donna Dea Holland . 2018. SAS® Certification Prep
#'             Guide: Statistical Business Analysis Using SAS®9. Cary, NC: SAS Institute Inc.
#'
#' @examples
#' values <- rnorm(100)
#' x <- skewCoeff(data = values,sd = sd(values))
skewCoeff <- function(data, sd) {

  avg <- mean(data)
  n <- length(data)
  skew <- sum(sapply(data, function(data){ (((data - avg) / sd) ^ 3)  })) * (n / ((n - 1) * (n - 2)))

  return(skew)
}

####################################################################################################
#' Adjusted Fisher-Pearson Excess Sample Kurtosis
#'
#' @param data The data for which kurtosis is computed (vector)
#' @param sd The population standard deviation, used to compute kurtosis (scalar)
#'
#' @return The kurtosis value (scalar)
#' @export
#'
#' @references Shreve, Joni N. and Donna Dea Holland . 2018. SAS® Certification Prep
#'             Guide: Statistical Business Analysis Using SAS®9. Cary, NC: SAS Institute Inc.
#'
#' @examples
#' values <- rnorm(100)
#' x <- kurtCoeff(data = values, sd = sd(values))
kurtCoeff <- function(data, sd) {

  avg <- mean(data)
  n <- length(data)
  sumZ <- sum(sapply(data, function(data){ ((data - avg) / sd) ^ 4 }))
  kurt <- ((n * (n + 1)) / ( (n - 1) * (n - 2) * (n - 3)  ) * sumZ) - ((3 * ((n - 1) ^ 2)) / ((n - 2) * (n - 3)))
  return(kurt)
}
#MISC################################################################################################
#' Converts Sample Standard Deviations into Population Equivalents
#'
#' This function converts a sample standard deviation (SD) input into the population equivalent. This code
#' is vectorized to convert several sample standard deviations for univariate distributions of identical
#' sample-sizes, if desired.
#'
#' @param s The sample SD(s) (vector)
#' @param n The sample-size for each SD to be converted (vector)
#'
#' @return The population SD(s) (vector)
#' @export
#'
#' @references Ruscio, J. (2021). Fundamentals of research design and statistical analysis. Ewing, NJ:
#'             The College of New Jersey, Psychology Department.
#' @examples
#' values <- rnorm(100)
#' x <- popSD(s = sd(values),n = 100)
popSD <- function(s, n) {

  sigma <- sapply(s, function(s) {  sqrt( ( (s ^ 2) * (n - 1) ) / n ) } )

  return(sigma)
}
