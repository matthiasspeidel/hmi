#' The function for hierarchical imputation of semicontinuous variables.
#'
#' The function is called by the wrapper. We consider data to be "semicontinuous" when
#' more than 5\% of the (non categorical) observations.\cr
#' For example in surveys a certain portion of people, when asked for their income,
#' report "0", which clearly violates the assumption of income to be (log-) normally distributed.
#' @param y_imp A Vector with the variable to impute.
#' @param X_imp A data.frame with the fixed effects variables.
#' @param heap A numeric value saying to which value the data might be heaped.
#' @param pvalue A numeric between 0 and 1 denoting the threshold of p-values a variable in the imputation
#' model should not exceed. If they do, they are excluded from the imputation model.
#' @param rounding_degrees A numeric vector with the presumed rounding degrees.
#' @return A n x 1 data.frame with the original and imputed values.
#' @export
imp_semicont_single <- function(y_imp,
                        X_imp,
                        heap = NULL,
                        pvalue = 0.2,
                        rounding_degrees = c(1, 10, 100, 1000)){

  if(is.null(heap)){
    heap <-  list_of_spikes_maker(data.frame(y_imp))$y_imp
  }

  if(is.null(heap)){
    heap <- 0
  }

  tmp_data <- cbind(y_imp, X_imp)
  n <- nrow(tmp_data)

  # transform y_imp into a binary variable,
  # with 0 representing a heaped value and 1 a non heaped value
  # NA values will remain NA.
  y_binary <- factor(rep(NA, length(y_imp)), levels = c(0, 1))

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
                                   X_imp = X_imp,
                                   pvalue = pvalue,
                                   rounding_degrees = rounding_degrees)


  # the data points where the binary varriable is "1" (meaning continuous)
  # are used for the continuous imputation
  y1_imp <- imp_cont_single(y_imp = y_imp[what_method == 1],
                            X_imp = X_imp[what_method == 1, , drop = FALSE],
                            pvalue = pvalue,
                            rounding_degrees = rounding_degrees)

  # set the final value of y:
  # the observations with method 1 (continuous (non hepead) observation)
  # get the continuously imputed values
  y_tmp <- array(NA, dim = length(y_imp))
  y_tmp[what_method == 1] <- y1_imp[, 1]

  # the observations with method 0 (heaped observation)
  # get the heap
  y_tmp[what_method == 0] <- heap

  y_ret <- data.frame(y_ret = y_tmp)

  return(y_ret)
}

