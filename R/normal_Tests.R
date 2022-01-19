###################################################################################
#' Kolmogorov-Smirnov-Lilliefors Test
#'
#' This function computes the Lilliefors variant of the one-sample Kolmogorov-Smirnov test.
#'
#' Molin & Abdi's (1998) algorithmic approximation of p-values is used for hypothesis-testing.
#' Note that this algorithm requires the imputation of 0.0 for negative output when p-values
#' would otherwise be low in value (< 0.001) using other methods. A similar issue with extremely
#' large values requires the imputation of 1.0 for values larger than 1.0 when p > .99.
#'
#' @param data The data of a univariate distribution for which the test statistic is computed (vector)
#' @param alpha The two-sided decision threshold used for hypothesis-testing (scalar)
#' @param j The # hypotheses tested; used to compute a Bonferonni correction, if applicable;
#'          should remain at its default if multiple testing is not an issue (scalar)
#' @param warn Used for printing a warning message when negative values are imputed to 0.0 (boolean)
#'
#' @return An object including the test statistic, p-value, and a significance flag (list)
#' @export
#'
#' @references Lilliefors, H.W. (1967). On the Kolmogorov-Smirnov Test for Normality with Mean and Variance Unknown. Journal of the American Statistical Association, 62, 399-402.
#'
#'             Molin, P., & Abdi, H. (1998). New Tables and numerical approximation for the KolmogorovSmirnov/Lillierfors/Van Soest test of normality.
#'
#' @examples
#' values <- rnorm(100)
#' x <- KSLTest(data = values)
KSLTest <- function(data, alpha = .05, j = 1, warn = T) {
  N <- length(data)
  alpha <- alpha / j

  phiCDF <- sort(pnorm(scale(data)))
  empCDF <- sort(ecdf(scale(data))(scale(data)))
  d0 <- abs(phiCDF - empCDF)
  dPrev <- abs(phiCDF - c(0,empCDF[-N]))
  d <- max(c(d0,dPrev))
  A <- (-( 1.30748185078790 + N) + sqrt( (1.30748185078790 + N)^2 - (4 *  0.08861783849346) *
		                                (0.37872256037043 - (1/d)^2)
                                       )
        ) / (2 *  0.08861783849346)
  AProb <- (-0.37782822932809 + (1.67819837908004 * A)
          + (-3.02959249450445 * A^2) + (2.80015798142101 * A^3)
          + (-1.39874347510845 * A^4) + (0.40466213484419 * A^5)
          + (-0.06353440854207 * A^6) + (0.00287462087623 * A^7)
          + (0.00069650013110 * A^8)  + (-0.00011872227037 * A^9)
          + (0.00000575586834 * A^10) )
 # AProb yields negative values when p is extremely small, requiring a floor:
 if (AProb < 0.0) {
   AProb <- 0.00
   if (warn)
     warning("Negative p-value imputed to 0.0")
 } else if (AProb > 1.0) {
   AProb <- 1.0
   if (warn)
     warning("P-value greater than 1.0 imputed to 1.0")
 }

 sigFL <- AProb < alpha
 names(sigFL) <- "Sig."
 names(d) <- "D"
 names(AProb) <- "P-value"
 stats <- list(round(d,3),round(AProb,3),sigFL)

 return(stats)
}
###################################################################################
#' Shapiro-Wilk Test
#'
#' This function is a wrapper for shapiro.test() from the stats package. Options
#' added include an ability to toggle a Bonferonni correction for significance, a
#' corresponding significance flag, and reorganized output to facilitate integration
#' with the Rita package.
#'
#' Note that when the sample-size of the input vector is > 5000, resampling with
#' replacement is used to proceed with hypothesis-testing with a vector of 5000 elements. 
#' When N < 3, testing is terminated.
#'
#' @param data Data of a univariate distribution for which the test statistic is computed
#'             (vector)
#' @param alpha The two-sided decision threshold used for hypothesis-testing
#' @param j The # hypotheses tested; used to compute a Bonferonni correction, if applicable;
#'          should remain at its default if multiple testing is not an issue (scalar)
#' @param warn Used for printing a warning message when resampling is performed on 
#'             sample-sizes > 5000 or when testing is terminated for N < 3 (boolean)
#'
#' @return An object including the test statistic, p-value, and a significance flag (list)
#' @export
#'
#' @references Patrick Royston (1982). An extension of Shapiro and Wilk's W test for normality to large samples. Applied Statistics, 31, 115--124. 10.2307/2347973
#'
#'             Patrick Royston (1982). Algorithm AS 181: The W test for Normality. Applied Statistics, 31, 176--180. 10.2307/2347986
#'
#'             Patrick Royston (1995). Remark AS R94: A remark on Algorithm AS 181: The W test for normality. Applied Statistics, 44, 547--551. 10.2307/2986146
#'
#' @examples
#' values <- rnorm(100)
#' x <- SWTest(data = values)
SWTest <- function(data, alpha = .05, j = 1, warn = T) {
  n <- length(data)
  if ( (n >= 3) & (n <= 5000)) {
    alpha <- alpha / j
    data <- sort(data)
    results <- shapiro.test(data)
    sigFL <- results[[2]] < alpha
    names(sigFL) <- "Sig."
    stats <- list(round(results[[1]],3),round(results[[2]],3),sigFL)
    names(stats[[2]]) <- "P-value"

    return(stats)
  } else if (n > 5000) {
    shuffled <- sample(data, size = 5000, replace = T)
    alpha <- alpha / j
    shuffled <- sort(shuffled)
    results <- shapiro.test(shuffled)
    sigFL <- results[[2]] < alpha
    names(sigFL) <- "Sig."
    stats <- list(round(results[[1]],3),round(results[[2]],3),sigFL)
    names(stats[[2]]) <- "P-value"
    
    if (warn) {
      warning("N > 5000: Resampling wih replacement for 5000 elements was performed 
               on the input data.")
    }

    return(stats)
  } else {
    if (warn)
      warning("Note: SW Test terminated due to N < 3.")

    return(list( rep("F",3) ))
  }
}
###################################################################################
#' Anderson-Darling Test
#'
#' This function computes the one-sample Anderson-Darling test statistic and p-value
#' for fit to a normal distribution.
#'
#' An adjusted statistic provided by D'agostino & Stephens (1986) is used,
#' where the mean and variance of the population are treated as unknown. D'agostino & Stephen's (1986)
#' text provides the equations used to obtain the function's p-values.
#'
#' @param data Data of a univariate distribution for which the test statistic is computed
#'            (vector)
#' @param alpha The two-sided decision threshold used for hypothesis-testing
#' @param j The # hypotheses tested; used to compute a Bonferonni correction, if applicable;
#'          should remain at its default if multiple testing is not an issue (scalar)
#'
#' @return An object including the test statistic, p-value, and a significance flag (list)
#' @export
#'
#' @references D'agostino, R. B., & Stephens, M. A. (1986). Goodness-of-fit-techniques (Vol. 68). CRC press.
#'
#' @examples
#' values <- rnorm(100)
#' x <- ADTest(data = values)
ADTest <- function(data, alpha = .05, j = 1) {
  alpha <- alpha / j
  data <- sort(data)
  N <- length(data)

  phiCDF <- sort(pnorm(scale(data)))
  logComb <- (sapply(1:N, function(x) { (2 * x) - 1 })* log(phiCDF)) +
			 (log(1 - phiCDF) * sapply(1:N, function(x) {2 * (N - x) + 1}))
  S <- sum(logComb) / N
  AD <- ((-N) - S) * (1 + (.75 / N) + (2.25 / N^2))
  if (AD <= .34) {
    p <- 1 - exp(-13.436 + (101.14 * AD) - (223.73 * AD^2))
  } else if ((AD > .20) & (AD < .34)) {
    p <- 1 - exp(-8.318 + (42.796 * AD) - (59.938 * AD^2))
  } else if ((AD > .34) & (AD < .60)) {
    p <- exp(0.9177 + -(4.279 * AD) + -(1.38 * AD^2))
  } else {
    p <- exp(1.2937 + -(5.709 * AD) + (0.0186 * AD^2))
  }

  sigFL <- p < alpha
  names(sigFL) <- "Sig."
  names(AD) <- "AD"
  names(p) <- "P-value"
  stats <- list(round(AD,3),round(p,3),sigFL)

  return(stats)
}
###################################################################################
#' D'agostino Pearson Omnibus Test
#'
#' This function computes the D'agostino Pearson omnibus test using adjusted Fisher-
#' Pearson skewness and kurtosis estimators.
#'
#' @param data Data of a univariate distribution for which the test statistic is computed
#'            (vector)
#' @param alpha The two-sided decision threshold used for hypothesis-testing
#' @param j The # hypotheses tested; used to compute a Bonferonni correction, if applicable;
#'          should remain at its default if multiple testing is not an issue (scalar)
#' @param warn Used for printing a warning message when testing is terminated for N < 8 (boolean)
#'
#' @return An object including the test statistic, p-value, and a significance flag (list)
#' @export
#'
#' @references D'agostino, R. B., & Stephens, M. A. (1986). Goodness-of-fit-techniques (Vol. 68). CRC press.
#'
#'             D’agostino, R. B., & Belanger, A. (1990). A Suggestion for Using Powerful and Informative Tests of Normality. The American Statistician, 44(4), 316–321. https://doi.org/10.2307/2684359
#'
#'             Shreve, Joni N. and Donna Dea Holland . 2018. SAS® Certification Prep
#'             Guide: Statistical Business Analysis Using SAS®9. Cary, NC: SAS Institute Inc.
#'
#' @examples
#' values <- rnorm(100)
#' x <- DPTest(data = values)
DPTest <- function(data, alpha = .05, j = 1, warn = T) {
  alpha <- alpha / j
  data <- sort(data)
  N <- length(data)
  popSD <- popSD(sd(data),N)

  # Compute the D'agostino test statistic for skew:
  c <- ( 3 * (N^2 + 27 * N -70) * (N + 1) * (N + 3) ) / ( (N - 2) * (N + 5) * (N + 7) * (N + 9) )
  wSq <- sqrt(2 * (c - 1)) - 1
  a <- sqrt((wSq - 1) / 2)
  b <- 1 / sqrt(log(sqrt(wSq)))
  mu <- a * skewCoeff(data,popSD) * sqrt( ((N + 1) * (N + 3)) / (6 * (N - 2)) )
  zSkew <- b * log(mu + sqrt(mu^2 + 1))

  # Compute the D'agostino test statistic for kurtosis:
  d <- sqrt( ((N + 1)^2 * (N + 3) * (N + 5)) / (24 * N * (N - 2) * (N - 3)) )
  e <- ((6 * (N^2 - 5 * N + 2)) / ((N + 7) * (N + 9))) * sqrt( (6 * (N + 3) * (N + 5)) / (N * (N - 2) * (N - 3)) )
  f <- (6 + (8 / e)) * (2 / e + sqrt(1 + 4 / e^2))
  g <- d * ((kurtCoeff(data,popSD) + 3) - ((3 * (N - 1)) / (N + 1)) ) * sqrt(2 / (f - 4))
  v <- (1 - (2 / f)) / (1 + g)
  r <- 2 / (9 * f)
  zKurt <- (1 - r - v^(1/3)) / sqrt(r)

  # Compute omnibus statistic:
  K <- zSkew^2 + zKurt^2

  # Compute p-value:
  p <- pchisq(K,df = 2, lower.tail = F)
  sigFL <- p < alpha
  names(sigFL) <- "Sig."
  names(K) <- "K"
  names(p) <- "P-value"
  stats <- list(round(K,3),round(p,3),sigFL)

  if (N >= 8) {
    return(stats)
  } else {
    if (warn)
      warning("Note: DP Test terminated due to N < 8.")

    return(list( rep("F",3) ))
  }
}
###################################################################################
#' Jarque-Bera Test
#'
#' This function performs the Jarque-Bera test for normality using adjusted Fisher-
#' Pearson skewness and kurtosis coefficients.
#'
#' Large samples (N >= 2000) use p-values obtained with reference to the chi-square
#' distribution, whereas smaller samples output p-values obtained via bootstrapping.
#' When N < 4, testing is terminated.
#'
#' @param data Data of a univariate distribution for which the test statistic is computed
#'            (vector)
#' @param alpha The two-sided decision threshold used for hypothesis-testing
#' @param j The # hypotheses tested; used to compute a Bonferonni correction, if applicable;
#'          should remain at its default if multiple testing is not an issue (scalar)
#' @param N_Sample The # samples used to generate the bootstrapped sampling distribution,
#'                 in cases when N < 2000 (scalar)
#' @param warn Used for printing a warning message when boostrapping is performed for 
#'             sample-sizes < 2000 or when testing is terminated for N < 4 (boolean)
#'
#' @return An object including the test statistic, p-value, and a significance flag (list)
#' @export
#'
#' @references Jarque, C. M. and Bera, A. K. (1980). Efficient test for normality, homoscedasticity and serial independence of residuals. Economic Letters, 6(3), pp. 255-259.
#'
#'             Shreve, Joni N. and Donna Dea Holland . 2018. SAS® Certification Prep
#'             Guide: Statistical Business Analysis Using SAS®9. Cary, NC: SAS Institute Inc.
#' @examples
#' values <- rnorm(100)
#' x <- JBTest(data = values)
JBTest <- function(data, alpha = .05, j =1, N_Sample = 10000, warn = T) {
  alpha <- alpha / j
  data <- sort(data)
  N <- length(data)
  popSD <- popSD(sd(data),N)

  # Compute p-value via chisq if N >= 2000:
  if (N >= 2000) {
    # Compute the test statistic:
    JB <- (N /6) * skewCoeff(data,popSD)^2 + (kurtCoeff(data,popSD)^2 / 4)

    # Compute p:
    p <- pchisq(JB,df = 2, lower.tail = F)

    sigFL <- p < alpha
    names(sigFL) <- "Sig."
    names(JB) <- "JB"
    names(p) <- "P-value"
    stats <- list(round(JB,3),round(p,3),sigFL)

    return(stats)
  # Compute p-value via simulation otherwise:
  } else if (N >= 4) {
    if (warn) {
      warning("N < 2000: Output p-values obtained via bootstrapping.")
    }

    # Compute the test statistic:
    JB <- (N /6) * skewCoeff(data,popSD)^2 + (kurtCoeff(data,popSD)^2 / 4)

    empDist <- rep(NA,N_Sample)
    parent <- rep(NA,N)
    for (i in 1:N_Sample) {
      parent <- data[sample(1:N,N, replace = T)]
      empDist[i] <- (N / 6) * skewCoeff(parent,popSD(sd(parent),N))^2 + (kurtCoeff(parent,popSD(sd(parent),N))^2 / 4)

    }

    # Compute p:
    p <- ecdf(empDist)(JB)

    sigFL <- p < alpha
    names(sigFL) <- "Sig."
    names(JB) <- "JB"
    names(p) <- "P-value"
    stats <- list(round(JB,3),round(p,3),sigFL)

    return(stats)
  #Terminate testing when N < 4:
  } else {
    if (warn)
      warning("Note: JB Test terminated due to N < 4.")

    return(list( rep("F",3) ))    
  }
}
###################################################################################
#' Chi-Square Test
#'
#' This function computes the chi-square test for normality.
#'
#' Bins are created by cutting the data to ensure that values within these intervals would be
#' equally probable if data are normal (Moore, 1986). By default, this function assumes
#' that all relevant parameters (mu, sigma) are estimators, fixing the degrees of freedom
#' at df = 3.
#'
#' @param data Data of a univariate distribution for which the test statistic is computed
#'            (vector)
#' @param alpha The two-sided decision threshold used for hypothesis-testing
#' @param j The # hypotheses tested; used to compute a Bonferonni correction, if applicable;
#'          should remain at its default if multiple testing is not an issue (scalar)
#' @param df The degrees of freedom used to test for significance against the sampling
#'           distribution (scalar)
#'
#' @return An object including the test statistic, p-value, and a significance flag (list)
#' @export
#'
#' @references Moore, D.S., (1986) Tests of the chi-squared type. In: D'agostino, R.B. and Stephens, M.A., eds.: Goodness-of-Fit Techniques. Marcel Dekker, New York.
#'
#' @examples
#' values <- rnorm(100)
#' x <- chisqTest(data = values)
chisqTest <- function(data, alpha = .05, j = 1, df = 3) {
  alpha <- alpha / j
  data <- sort(data)
  N <- length(data)
  mu <- mean(data)
  sigma <- sd(data)

  bins <- cut(data,breaks = ceiling(2 * N^(2/5)), include.lowest = T)
  binsX <- trimws(chartr(old = "[",new = " ",paste0(levels(bins), collapse = " ")))
  binsX <- chartr(old = "(", new = " ",binsX)
  binsX <- chartr(old = "]", new = " ",binsX)
  binsX <- chartr(old = ",", new = " ",binsX)
  binsX <- strsplit(binsX,split = "  ")
  binsX <- unlist(strsplit(unlist(binsX),split = " "))
  binsX <- as.numeric(binsX[binsX != ""])
  lowerLim <- binsX[seq(1,length(binsX),by = 2)]
  upperLim <- binsX[seq(2,length(binsX),by = 2)]


  E <- rep(N / ceiling(2 * N^(2/5)),length(lowerLim))
  chiSq <- 0
  for (i in 1:length(lowerLim)) {
    if (i != 1) {
      O <- eval(data > lowerLim[i])
      O <- sum(eval(data[O] <= upperLim[i]))
    } else if (i == 1) {
      O <- eval(data >= lowerLim[i])
      O <- sum(eval(data[O] <= upperLim[i]))
    }
    chiSq <- chiSq + ((O - E[i])^2 / E[i])
  }
  p <- pchisq(chiSq,df = 3, lower.tail = F)

  sigFL <- p < alpha
  names(sigFL) <- "Sig."
  names(chiSq) <- "chiSq"
  names(p) <- "P-value"
  stats <- list(round(chiSq,3),round(p,3),sigFL)

  return(stats)
}
