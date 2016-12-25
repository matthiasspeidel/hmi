#' The function for hierarchical imputation of semicontinous variables.
#'
#' The function is called by the wrapper. We consider data to be "semicontinuous" when
#' more than 5\% of the (non categorical) observations.\cr
#' For example in surveys a certain portion of people, when asked for their income,
#' report "0", which clearly violates the assumption of income to be (log-) normally distributed.
#' @param y_imp_multi A Vector with the variable to impute.
#' @param X_imp_multi A data.frame with the fixed effects variables.
#' @param heap A scalar saying to which value the data might be heaped.
#' @param M An integer defining the number of imputations that should be made.
#' @return A n x M matrix. Each column is one of M imputed y-variables.
imp_semicont<- function(y_imp_multi,
                        X_imp_multi,
                        heap = 0,
                        M = 10){


  tmp_data <- cbind(y_imp_multi, X_imp_multi)
  n <- nrow(tmp_data)

  #the missing indactor indicates, which values of y are missing.
  mis_indicator <- is.na(y_imp_multi)
  #get the defaults values for heap
  if(is.null(heap)) heap <- 0

  # transform y_imp_multi into a binary variable,
  # with 0 representing a heaped value and 1 a non heaped value
  # NA values will remain NA.
  y_binary <- y_imp_multi

  # The observations beeing heaped and not NA...
  condition_0 <- (y_imp_multi == heap) & !is.na(y_imp_multi)
  #... are set to be 0
  y_binary[condition_0] <- 0

  # The observations beeing not heaped and not NA...
  condition1 <- (y_imp_multi != heap) & !is.na(y_imp_multi)
  #... are set to be 1
  y_binary[condition1] <- 1

  # Use the imputation function of the binary variable on the indicator
  # to figure out whether a missing value shall get the value of the heap or
  # a continuous.
  # For the data points with an observed y_imp_multi,
  # this also indicates whether they are used in the continuous imputation model or not.
  what_method <- imp_binary(y_imp_multi = y_binary,
                            X_imp_multi = X_imp_multi)

  # set up the results matrix
  y_imp <- array(NA, dim = c(n, M))
  for(i in 1:M){

    # the data points where the binary varriable is "1" (meaning continuous)
    # are used for the continuous imputation
    y1_imp <- imp_cont(y_imp_multi = y_imp_multi[what_method[, i] == 1],
                              X_imp_multi = X_imp_multi[what_method[, i] == 1, , drop = FALSE])

    # set up the final i-th imputation vector
    y_tmp <- what_method[, i]

    # the data points where the binary imputation said, that they shall be (or already are) continuous,
    # get the value of the continuous imputation
    y_tmp[what_method[, i] == 1] <- y1_imp

    y_imp[ , i] <- y_tmp

  }

  return(y_imp)
}

