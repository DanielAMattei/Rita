#' Rita
#'
#' R Exploratory Data Analysis (REDA; pronounced "rita")
#' summarizes an input dataset by the M, SD + 5-number summary + third and fourth moments
#' and visualizes the data according to an algorithm or as specified by the user.
#' In addition, Rita will provide the results of one or several normality tests.
#' Lastly, Rita normalizes the dataset with several methods and provides
#' visualizations of the best performing method to the user.
#'
#' Any rows with missing values (NAs) are removed for calculation purposes; if desired,
#' incomplete records should be imputed prior to calling Rita. In addition, note that any
#' columns not numeric type or coercible to numeric are by excluded from analysis.
#'
#' @param data Input dataset (matrix, dataframe, or vector).
#'             For a univariate distribution, submit a vector or a subsetted
#'             matrix or dataframe. If results for many univariate distributions are
#'             desired, submit a matrix or dataframe with each column representing a
#'             given variable if all distributions are of the same sample-size. If not,
#'             it is recommended to call Rita repeatedly for each variable.
#' @param test Desired normality test (scalar). By default (test = 1), Rita will
#'             present the results of the Shapiro-wilk test to the user.
#'
#'             test = 1: Shapiro-Wilk (SW)
#'
#'             test = 2: Kolmogorov-Smirnov/Lilliefors (KSL)
#'
#'             test = 3: Anderson-Darling (AD)
#'
#'             test = 4: Jarque-Bera (JB)
#'
#'             test = 5: D'Agostino Pearson Omnibus (DP)
#'
#'             test = 6: Chi-square test (chiSq)
#'
#'             test = 7: Results of all tests for the best performing transformation
#'
#'             The order of the tests printed corresponds to the order of the variables stored within
#'             the input dataset.
#' @param xform Desired normalization method (scalar). By default (xform = 1), Rita
#'              will assess which method performs best and (a.) return the transformed data to
#'              the user, and (b.) visualize the data according to the settings of the plot argument.
#'
#'              Please note that, per the recommendations of Osborne (2002), a constant is added
#'              prior to logarithmic and inverse transformations to ensure that the minimum value
#'              is anchored at 1, and prior to the square-root transformation to ensure a left anchor
#'              of 0.
#'
#'              Similarly, the arc-sine and logit transformations are applied after converting the units,
#'              if needed, to ensure that variables are bounded between 0 and 1.
#'
#'      	      The "best performing" method is identified by comparing goodness-of-fit to the
#'              straight line of the QQ plot for the quantiles of the data normalized by a given method
#'              and the standard normal distribution. If a tie is present between transformations for a
#'              variable, one of the best performing transformations is arbitrarily selected.
#'
#'      	      xform = 1: Best performing method is presented (excluding the Rankit)
#'
#'          	  xform = 2: Logarithmic transform
#'
#'          	  xform = 3: Inverse/reciprocal transform
#'
#'          	  xform = 4: Square-root transform
#'
#'          	  xform = 5: Arc-sine transform
#'
#'              xform = 6: Logit transform
#'
#'          	  xform = 7: Rankit transform
#' @param alpha The two-sided decision threshold used for normality hypothesis-testing (scalar)
#' @param j The # hypotheses tested; used to compute a Bonferonni correction, if applicable;
#'          should remain at its default if multiple testing is not an issue (scalar)
#' @param autoPlot Desired plotting method (boolean). By default (plot = 1),
#'                 the visualization will be implicitly chosen based on
#'                 extracted features of the dataset.
#'
#'                 When autoPlot = F, values of additional plotting arguments are used to
#'                 determine the visualizations provided to the user.
#'
#'            	   When autoPlot = T:
#'
#'	               Histograms are always generated for discrete data.
#'
#'	               Density plots are always generated for continuous data.
#'
#'                 Strip plots are generated when the # distinct values
#'                 are <= 20 AND the # datapoints are 15 <= x <= 150.
#'
#'	               Violin plots are instead generated in lieu of the strip plots
#'                 created when the above conditions are not met.
#'
#'                 Lastly, density plots for each (transformed*) variable are generated.
#'
#'                 *Transformed variables correspond to the choice made by the user for the xform argument
#'                 or to the best-performing transformation for each variable when xform = 1.
#'
#'                 All plots are drawn in the R console and saved as plotting objects.
#'
#' @param histPlot Whether to generate histograms for each variable (boolean).
#' @param densPlot Whether to generate density plots for each variable (boolean).
#' @param stripPlot Whether to draw strip plots for each variable (boolean).
#' @param violinPlot Whether to draw violin plots for each variable (boolean).
#' @param xformPlot Whether to draw density plots for each transformed variable (boolean).
#' @param seed Number used for reproduction of random number generator results (scalar)
#'
#' @return An object containing the dataset of the best performing transformation for
#'         each variable and the specified plots (list)
#' @export
#'
#' @examples
#' values <- rnorm(100)
#' x <- Rita(data = values)
Rita <- function(data, test = 1, xform = 1, alpha = .05, j =1, autoPlot = T, histPlot = F, densPlot = F,
                 stripPlot = F, violinPlot = F, xformPlot = F, seed = 10) {
  # Set a seed:
  set.seed(seed)

  # Coerce data to dataframes to generalize code for NA-checking:
  postData <- as.data.frame(data)

  # Store # NAs in each column + row #s with incomplete records:
  NAs <- apply(postData,2,function(postData) { sum(eval(is.na(postData) == T)) })
  NARows <- nrow(postData[complete.cases(postData) == F,])
  if (is.null(NARows))
    NARows <- 0

  # Store values of columns coerced to numeric to distinguish non-numeric data:
  postData <- suppressWarnings(apply(postData,2,as.numeric))

  # Omit entire columns composed of NAs  + remove cols with uniform distributions:
  NAFL <- uniFL <- rep(0,ncol(postData))
  for (i in 1:ncol(postData)) {
    colFL <- sum(eval(is.na(postData[,i] == T)))
    if (colFL == length(postData[,i]))
      NAFL[i] <- i
    if ((length(unique(postData[,i])) == 1) & (colFL != length(postData[,i])))
      uniFL[i] <- i
  }
  if (sum(NAFL) > 0)
    postData <- postData[,-NAFL]
  if (sum(uniFL) > 0) {
    initVarNames <- colnames(postData)[uniFL]
    postData <- postData[,-uniFL]
  }
  # Omit (at least) partially incomplete records:
  postData <- as.data.frame(na.omit(postData))
  univNRow <- nrow(postData)
  univNCol <- ncol(postData)

  SD <- apply(postData,2,sd)
  # Obtain the population SD to compute 3rd + 4th moments:
  popSD <- popSD(SD,univNRow)
  skew <- mapply(skewCoeff,as.list(postData),popSD)
  kurt <- mapply(kurtCoeff,as.list(postData),popSD)

  # Create list object of descriptive statistics for the input data:
  preStats <- round(rbind(as.data.frame(apply(postData,2,summary)),SD,skew,kurt),3)
  dimnames(preStats)[[1]][7:9] <- c("SD","Skewness","Kurtosis")

  # Perform transformations based on the xform argument:
  xformAgg <- xformNames <- vector(length = univNCol, mode = "list")
  for (i in 1:univNCol) {
    xformAgg[[i]] <- MasterXform(xform,postData[,i])
  }
  xformAgg <- do.call(cbind,xformAgg)
  xformBest <- matrix(NA,nrow(xformAgg),univNCol)
  xformNames <- if (xform != 1) {
                  c("Log","Inverse","Square-root","Arc-sine","Logit","Rankit")[xform - 1]
                } else {
                  rep(c("Log","Inverse","Square-root","Arc-sine","Logit"),univNCol)
                }
  xformBestNames <- vector(mode = "list")

  # Calculate R^2 for all normalizations via the QQ plot:
  QQOutput <- apply(xformAgg,2, function(xformAgg) { qqnorm(xformAgg, plot.it = F) })
  QQEst <- rep(NA,length(QQOutput))

  # Correlations are stored such that transformations are listed sequentially with variables held
  # constant:
  QQEst <- sapply(QQOutput,function(QQOutput) { cor(QQOutput$x,QQOutput$y) })
  QQEst <- QQEst^2

  # Identify the best performing transformations by R^2:
  # Best R^2 for each variable is stored below:
  maxEst <- rep(NA,univNCol * 5)
  if (xform == 1) {
     lower <- seq(1,10e5, by = 5)
     upper <- seq(5,univNCol * 5, by = 5)
     for (i in 1:(univNCol * 5)) {
       maxEst[i] <- max(QQEst[lower[i]:upper[i]])
       # Subset transformed data with best R^2 for normality testing: + visuals
       xformBest[,i] <- xformAgg[,which(maxEst[i] == QQEst)[1]]
       xformBestNames[[i]] <- xformNames[which(maxEst[i] == QQEst)]
       if (upper[i] == max(upper))
         break
     }
  } else {
    maxEst <- QQEst
    xformBest <- xformAgg
    xformBestNames <- xformNames
  }
  maxEst <- na.omit(maxEst)

  # Conduct normality tests based on the test argument:
  norColNames <- c("SW","KSL","AD","JB","DP","chiSq")
  norRowNames <- c("Stat.","P-Value","Sig.")
  if (test != 7) {
    norOutput <- as.data.frame(matrix(,3,univNCol))
    for (i in 1:univNCol) {
      norOutput[,i] <- unlist(MasterTest(test,xformBest[,i],alpha,j))
	  norOutput[3,i] <- if (norOutput[3,i]) {
						  "T"
						} else {
						  "F"
						}
    }
    rownames(norOutput) <- norRowNames
    colnames(norOutput) <- rep(norColNames[test],univNCol)
    norAgg <- norOutput
  } else {
    norOutput <- as.data.frame(matrix(,3,(test - 1)))
    norAgg <- list(rep(NA,univNCol * (test - 1)))
    for (i in 1:univNCol) {
      for (j in 1:(test - 1)) {
        norOutput[,j] <- unlist(MasterTest(j,xformBest[,i],alpha,j))
	    norOutput[3,j] <- if (norOutput[3,i]) {
						    "T"
						  } else {
							"F"
						  }
      }
      rownames(norOutput) <- norRowNames
      colnames(norOutput) <- norColNames
      norAgg[[i]] <- norOutput
    }
  }

  # Prepare naming variables for final output objects:
  varNames <- names(postData)
  xformNames <- mapply(paste,varNames,xformBestNames)
  colnames(xformBest) <- xformNames
  names(norAgg) <- varNames[1:length(norAgg)]

  # Plot data according to the plot argument:
  if (autoPlot == T) {
    disDists <- conDists <- vector(length = univNCol, mode = "list")
    disInd <- conInd <- c()
    disFL <- 0
    conFL <- 0
    disLoc <- conLoc <- rep(NA,univNCol)
    suppPlot <- mainPlot <- vector(mode = "list")
    for (i in 1:univNCol) {
	  # Select for discrete data:
      if (all.equal(postData[,i],as.integer(postData[,i])) == T) {
		disDists[[i]] <- postData[,i]
        disFL <- 1
		disLoc[i] <- i
        disInd <- c(disInd,i)
      # Select for continuous data:
      } else {
        conDists[[i]] <- postData[,i]
		conFL <- 1
	    conLoc[i] <- i
        conInd <- c(conInd,i)
      }
      # Select for small # distinct values and small n +
      # Generate strip plots, if applicable:
      if ((length(unique(postData[,i])) <= 20) & (univNRow <= 150)) {
        suppPlot[[i]] <- stripplot(~postData[,i], xlab = varNames[i], jitter.data = T,
					             factor = 13)
      } else {
      # Generate violin plots, if applicable:
        suppPlot[[i]] <- bwplot(~postData[,i], xlab = varNames[i], panel = function(...) {
																		   panel.violin(...)
																		   panel.bwplot(...)})
      }
    }
	disLoc <- na.omit(disLoc)
	conLoc <- na.omit(conLoc)
    disDists <- disDists[disInd]
    conDists <- conDists[conInd]
    # Generate histograms, if applicable:
    if (disFL == 1) {
      hists <- vector(mode = "list")
      histNames <- varNames[disLoc]
      for (i in 1:length(disDists)) {
        hists[[i]] <- histogram(~disDists[[i]],xlab = histNames[i])
      }
      mainPlot <- c(mainPlot,hists)
    }
    # Generate density plots, if applicable:
    if (conFL == 1) {
      dens <- vector(mode = "list")
	  denNames <- varNames[conLoc]
      for (i in 1:length(conDists)) {
        dens[[i]] <- densityplot(~conDists[[i]],xlab = denNames[i])
      }
      mainPlot <- c(mainPlot,dens)
    }
    # Generate density plots for transformed variables:
    densXform <- vector(mode = "list")
    for (i in 1:univNCol) {
      densXform[[i]] <- densityplot(~xformBest[,i],xlab = xformNames[i])
    }
    plots <- list(c(rbind(mainPlot,densXform)),suppPlot)
  } else {
    hists <- dens <- strips <- violins <- densXform <- vector(mode = "list")
    for (i in 1:univNCol) {
      if (histPlot == T)
        hists[[i]] <- histogram(~postData[,i],xlab = varNames[i])
      if (densPlot == T)
        dens[[i]] <- densityplot(~postData[,i],xlab = varNames[i])
      if (stripPlot == T) {
        strips[[i]] <- stripplot(~postData[,i], xlab = varNames[i], jitter.data = T,
				                 factor = 13)
      }
      if (violinPlot == T) {
        violins[[i]] <- bwplot(~postData[,i], xlab = varNames[i], panel = function(...) {
  							                                                panel.violin(...)
																	        panel.bwplot(...)})
      }
      if (xformPlot == T)
        densXform[[i]] <- densityplot(~xformBest[,i],xlab = xformNames[i])
    }
    plots <- list(hists,dens,strips,violins,densXform)
  }

  # Display results:
  if (sum(histPlot,densPlot,stripPlot,violinPlot,xformPlot,autoPlot) != 0) {
    cat("\n")
    print(plots)
  } else {
    cat("\n")
    message("Note: Plots are turned off.")
  }
  if (sum(uniFL) > 0) {
    cat("\n")
    message("Note: Variables ",paste(initVarNames,collapse = ", ")," were removed due to uniform values.")
  }
  NAOutput <- paste("\nNo. missing values omitted from each column:", paste(dimnames(postData)[2][[1]],NAs, collapse = " ", sep = ": "),
					paste("\nNo. rows omitted with missing values:", NARows, "\n", collapse = " "), collapse = " ", sep = "\n")
  cat(NAOutput,"\n")
  cat("Desc. Stats of Raw Variable(s):\n")
  print(preStats)
  cat("\nNormality test results:\n")
  print(norAgg)
  if (test != 7)
  cat("Test type:", norColNames[test],"\n\n")

  return(list(xformBest,plots))
}

####################################################################################
#' Master Transformation Function
#'
#' This is a master function used to perform the appropriate transformation(s) within
#' the 'Rita' function.
#'
#' @param c Input specifying the test to run (scalar)
#' @param data The data of a univariate distribution for which the test statistic is computed (vector)
#'
#' @return Output from the appropriate subfunction (list)
#' @export
#'
#' @examples
#' values <- rnorm(100)
#' x <- MasterXform(c = 2, data = values)
MasterXform <- function(c, data) {
  switch(c,
    mapply(MasterXform,2:6,list(data)),
    logXform(data),
    inverseXform(data),
    squareXform(data),
    arcsineXform(data),
    logitXform(data),
    rankitXform(data)
  )
}
####################################################################################
#' Master Normality Testing Function
#'
#' This is a master function to call the appropriate test(s) to be used in the 'Rita' function.
#'
#' @param c Input specifying the test to run (scalar)
#' @param data The data of a univariate distribution for which the test statistic is computed (vector)
#' @param alpha The two-sided decision threshold used for hypothesis-testing (scalar)
#' @param j The # hypotheses tested; used to compute a Bonferonni correction, if applicable;
#'          should remain at its default if multiple testing is not an issue (scalar)
#'
#' @return An results object specific to the test designated with the 'c' argument (list)
#' @export
#'
#' @examples
#' values <- rnorm(100)
#' x <- MasterTest(c = 1, data = values)
MasterTest <- function(c, data, alpha = .05, j = 1) {
  switch(c,
    SWTest(data, alpha, j),
    KSLTest(data, alpha, j),
    ADTest(data, alpha, j),
    JBTest(data, alpha, j),
    DPTest(data, alpha, j),
    chisqTest(data, alpha, j)
    )
}
