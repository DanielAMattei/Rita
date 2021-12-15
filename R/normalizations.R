####################################################################################
#' Rankit Transformation
#'
#' This function transforms data via the Rankit, a member of the families of 'rank-based normalization
#' methods' and 'empirical normal quantile transformations' employed in both the social sciences
#' and quantitative genetics.
#'
#' @param sample The input data (vector)
#'
#' @return The Rankit-transformed data (vector)
#' @export
#'
#' @references Soloman, S. R., & Sawilowsky, S. S. (2009). Impact of rank-based normalizing transformations on the accuracy of test scores. Journal of Modern Applied Statistical Methods, 8(2), 9.
#'
#'             Peng, B., Robert, K. Y., DeHoff, K. L., & Amos, C. I. (2007, December). Normalizing a large number of quantitative traits using empirical normal quantile transformation. In BMC proceedings (Vol. 1, No. 1, p. S156). BioMed Central. doi: 10.1186/1753-6561-1-s1-s156
#'
#'             Bliss, C. I., Greenwood, M. L., & White, E. S. (1956). A rankit analysis of paired comparisons for measuring the effect of sprays on flavor. Biometrics, 12(4), 381-403.
#' @examples
#' values <- rnorm(100)
#' x <- rankitXform(values)
rankitXform <- function(sample) {
  sample <- qnorm((rank(sample) - .5) / length(sample))
  return(sample)
}
####################################################################################
#' Logarithmic Transformation
#'
#' This function imputes minimum values per the recommendations of Osborne (2002) and
#' subsequently transforms the data to a base-10 logarithmic scale.
#'
#' @param sample The input data (vector)
#'
#' @return The log-transformed data (vector)
#' @export
#'
#' @references Osborne, J. W. (2002). Notes on the use of data transformations. Practical Assessment, Research and Evaluation, 9(1), 42-50.
#'
#'             Osborne, J. W. (2002). The Effects of Minimum Values on Data Transformations. Retrieved from https://files.eric.ed.gov/fulltext/ED463313.pdf
#'
#' @examples
#' values <- rnorm(100)
#' x <- logXform(values)
logXform<- function(sample) {
   if (min(sample) < 1) {
     sample <- sample + (1 - min(sample))
   } else if (min(sample) > 1) {
     sample <- sample - (min(sample) - 1)
   }
   sample <- log10(sample)
   return(sample)
}
####################################################################################
#' Inverse/Reciprocal Transformation
#'
#' This function imputes minimum values per the recommendations of Osborne (2002) and
#' subsequently transforms the data using the reciprocal.
#'
#' @param sample The input data (vector)
#'
#' @return The reciprocal-transformed data (vector)
#' @export
#'
#' @references Osborne, J. W. (2002). Notes on the use of data transformations. Practical Assessment, Research and Evaluation, 9(1), 42-50.
#'
#'             Osborne, J. W. (2002). The Effects of Minimum Values on Data Transformations. Retrieved from https://files.eric.ed.gov/fulltext/ED463313.pdf
#'
#'
#' @examples
#' values <- rnorm(100)
#' x <- inverseXform(values)
inverseXform <- function(sample) {
  sample <- sample * -1
   if (min(sample) < 1) {
     sample <- sample + (1 - min(sample))
   } else if (min(sample) > 1) {
     sample <- sample - (min(sample) - 1)
   }
  sample <- 1 / sample
  return(sample)
}
####################################################################################
#' Square-root Transformation
#'
#' This function left anchors the minimum value to 0 per the recommendations of Osborne (2002)
#' and subsequently transforms the data by taking the square-root of each value.
#'
#' @param sample The input data (vector)
#'
#' @return The square-transformed data (vector)
#' @export
#'
#' @references Osborne, J. W. (2002). Notes on the use of data transformations. Practical Assessment, Research and Evaluation, 9(1), 42-50.
#'
#'             Osborne, J. W. (2002). The Effects of Minimum Values on Data Transformations. Retrieved from https://files.eric.ed.gov/fulltext/ED463313.pdf
#'
#' @examples
#' values <- rnorm(100)
#' x <- squareXform(values)
squareXform <- function(sample) {
   if (min(sample) < 0) {
     sample <- sample + (0 - min(sample))
   } else if (min(sample) > 0) {
     sample <- sample - (min(sample) - 0)
   }
  sample <- sqrt(sample)
  return(sample)
}
####################################################################################
#' Arcsine Transformation
#'
#' This function transforms the scale, if needed, to values of unity.
#' Then, the data is transformed by taking the arcsine of each value.
#' Per the recommendations of Osborne(2002), data points are left-anchored
#' at 0 to maximize the efficacy of the square-root transformation used enroute
#' to the arcsine.
#'
#' @param sample The input data (vector)
#'
#' @return The arcsine-transformed data (vector)
#' @export
#'
#' @references Osborne, J. W. (2002). Notes on the use of data transformations. Practical Assessment, Research and Evaluation, 9(1), 42-50.
#'
#'             Osborne, J. W. (2002). The Effects of Minimum Values on Data Transformations. Retrieved from https://files.eric.ed.gov/fulltext/ED463313.pdf
#'
#' @examples
#' values <- rnorm(100)
#' x <- arcsineXform(values)
arcsineXform <- function(sample) {
  if (min(sample) < 0) {
     sample <- sample + (0 - min(sample))
  }
  if (max(sample) > 1) {
    sample <- sample / max(sample)
  }
  sample <- asin(squareXform(sample))
  return(sample)
}
####################################################################################
#' Logit/Log-Odds Transformation
#'
#' This function transforms data via the logit/log-odds transformation.
#'
#' Initially, features of the input data are extracted and used to determine
#' an initial transformation to perform.
#'
#' All forms of data representing an underlying discrete scale are converted to proportions
#' of the total sample size, if needed. In these cases, values should be stored such that
#' elements are in absolute frequency, relative frequency, or percentage form.
#'
#' For non-count data, variables are shifted and bounded at [0,1] in a manner analogous to
#' the potential transformations of the scale performed by arcsineXform() prior to the arcine,
#' although transformed values are not expected to outperform more suitable transformations.
#'
#' Then, the empirical logit transformation is applied to avoid zeroes or ones, and the data
#' are transformed by taking the log-odds/logit of each value.
#'
#' @param sample The input data (vector, matrix, or dataframe)
#'
#' @return The logit-transformed data (vector)
#' @export
#'
#' @references Stevens, S., Valderas, J. M., Doran, T., Perera, R., & Kontopantelis, E. (2016). Analysing indicators of performance, satisfaction, or safety using empirical logit transformation. bmj, 352.
#'
#'             Osborne, J. W. (2002). Notes on the use of data transformations. Practical Assessment, Research and Evaluation, 9(1), 42-50.
#'
#'             Osborne, J. W. (2002). The Effects of Minimum Values on Data Transformations. Retrieved from https://files.eric.ed.gov/fulltext/ED463313.pdf
#'
#' @examples
#' values <- rnorm(100)
#' x <- logitXform(values)
logitXform <- function(sample) {
 # Coerce tables stored in arrays to vector form:
 sample <- as.vector(sample)

 total <- sum(sample)

 # Detect proportions if logic fails below -> Leave intact:
 # Detect non-count data ->  Bound at [0,1]:
 if (min(sample) < 0) {
   sample <- sample + (0 - min(sample))
 }
 if (total != 1.0 & total != 100) {
   if (max(sample) >= 1)
     sample <- sample / max(sample)
 # Detect absolute frequencies or percentages -> Convert to proportions
 } else if (total == 100) {
   if ( sum(sapply(sample,function(sample) sample%%1) == 0) == length(sample)) {
     sample <- sample / sum(sample)
   } else {
     sample <- sample / 100
   }
 }
 # Replace zeroes and ones with a variant of the
 # empirical logit transformation, if needed:
 epsilon <- .5 * unique(sort(sample))[2]
 sample[sample == 0] <- .5 * unique(sort(sample))[2]
 sample[sample == 1] <- 1 - epsilon
 sample <- qlogis(sample)

 return(sample)
}
