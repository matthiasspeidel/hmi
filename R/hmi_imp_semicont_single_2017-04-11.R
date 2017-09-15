#' The function for hierarchical imputation of semicontinuous variables.
#'
#' The function is called by the wrapper. We consider data to be "semicontinuous" when
#' more than 5\% of the (non categorical) observations.\cr
#' For example in surveys a certain portion of people, when asked for their income,
#' report "0", which clearly violates the assumption of income to be (log-) normally distributed.
#' @param y_imp A Vector with the variable to impute.
#' @param X_imp A data.frame with the fixed effects variables.
#' @param heap A scalar saying to which value the data might be heaped.
#' @return A n x 1 data.frame with the original and imputed values.
imp_semicont_single <- function(y_imp,
                        X_imp,
                        heap = 0){


  tmp_data <- cbind(y_imp, X_imp)
  n <- nrow(tmp_data)

  #the missing indactor indicates, which values of y are missing.
  mis_indicator <- is.na(y_imp)
  #get the defaults values for heap
  if(is.null(heap)) heap <- 0

  # transform y_imp into a binary variable,
  # with 0 representing a heaped value and 1 a non heaped value
  # NA values will remain NA.
  y_binary <- y_imp

  # The observations beeing heaped and not NA...
  condition_0 <- (y_imp == heap) & !is.na(y_imp)
  #... are set to be 0.
  y_binary[condition_0] <- 0

  # The observations beeing not heaped and not NA...
  condition1 <- (y_imp != heap) & !is.na(y_imp)
  #... are set to be 1.
  y_binary[condition1] <- 1

  # Use the imputation function of the binary variable on the indicator
  # to figure out whether a missing value shall get the value of the heap or
  # a continuous.
  # For the data points with an observed y_imp,
  # this also indicates whether they are used in the continuous imputation model or not.
  what_method <- imp_binary_single(y_imp = y_binary,
                            X_imp = X_imp)


  # the data points where the binary varriable is "1" (meaning continuous)
  # are used for the continuous imputation
  y1_imp <- imp_cont_single(y_imp = y_imp[what_method == 1],
                            X_imp = X_imp[what_method == 1, , drop = FALSE])

  # set up the final i-th imputation vector
  y_tmp <- data.frame(what_method)

  # the data points where the binary imputation said, that they shall be (or already are) continuous,
  # get the value of the continuous imputation
  y_tmp[what_method == 1, 1] <- y1_imp

  y_ret <- data.frame(y_imp = y_tmp)

  return(y_ret)
}

