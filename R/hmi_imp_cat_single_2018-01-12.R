#' The function to impute unordered categorical variables
#'
#' The function uses regression trees for imputation implemented in \code{mice}.
#' The principle is the following:
#' For each observation it is calculated at which leave it would end.
#' Then one (randomly selected) observation of the other observations found on this leave
#' functions as a donor.
#' @param y_imp A Vector with the variable to impute.
#' @param X_imp A data.frame with the fixed effects variables.
#' @param rounding_degrees A numeric vector with the presumed rounding degrees.
#' @return A n x 1 data.frame with the original and imputed values as a factor.
imp_cat_single <- function(y_imp,
                         X_imp,
                         rounding_degrees = c(1, 10, 100, 1000)){

  if(min(table(y_imp)) < 2) {
    stop("Too few observations per category in a categorical target variable.")
  }

  # ----------------------------- preparing the X data ------------------
  # remove excessive variables
  X_imp <- cleanup(X_imp)

  # standardize X
  X_imp_stand <- stand(X_imp, rounding_degrees = rounding_degrees)


  #the missing indactor indicates, which values of y are missing.
  missind <- is.na(y_imp)

  n <- length(y_imp)

  #starting model
  ph <- sample_imp(y_imp)[, 1]
  tmp_0_all <- data.frame(target = ph)
  xnames_0 <- paste("X", 1:ncol(X_imp_stand), sep = "")
  tmp_0_all[xnames_0] <- X_imp_stand
  tmp_0_sub <- tmp_0_all[!missind, , drop = FALSE]

  reg_1_all <- nnet::multinom(target ~ 0 + ., data =  tmp_0_all, trace = FALSE)
  reg_1_sub <- nnet::multinom(target ~ 0 + ., data =  tmp_0_sub, trace = FALSE)

  X_model_matrix_1_all <- stats::model.matrix(reg_1_all)
  xnames_1 <- paste("X", 1:ncol(X_model_matrix_1_all), sep = "")


  #remove unneeded variables
  unneeded <- apply(stats::coefficients(reg_1_sub), 2, function(x) any(is.na(x)))
  xnames_2 <- xnames_1[!unneeded]

  tmp_2_all <- data.frame(target = as.factor(y_imp)) # mice needs the variable as a factor
  tmp_2_all[, xnames_2] <- X_model_matrix_1_all[, !unneeded, drop = FALSE]

  everything <- mice::mice(data = tmp_2_all, m = 1,
                     method = "cart",
                     predictorMatrix = (1 - diag(1, ncol(tmp_2_all))),
                     visitSequence = (1:ncol(tmp_2_all))[apply(is.na(tmp_2_all),2,any)],
                     post = vector("character", length = ncol(tmp_2_all)),
                     defaultMethod = "cart",
                     maxit = 10,
                     diagnostics = TRUE,
                     printFlag = FALSE,
                     seed = NA,
                     imputationMethod = NULL,
                     defaultImputationMethod = NULL,
                     data.init = NULL)

  #Initialising the returning vector
  y_ret <- y_imp

  y_ret[missind] <- everything$imp[[1]][, 1]

  return(data.frame(y_ret = factor(y_ret)))
}
