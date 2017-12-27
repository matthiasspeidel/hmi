#' The function for imputation of continuous variables.
#'
#' The function is called by the wrapper (hmi). It uses \code{mice} with the method "norm".
#' @param y_imp A Vector with the variable to impute.
#' @param X_imp A data.frame with the fixed effects variables.
#' @return A n x 1 data.frame with the original and imputed values.
imp_cont_single <- function(y_imp,
                      X_imp){


  #the missing indactor indicates, which values of y are missing.
  missind <- is.na(y_imp)
  # ----------------------------- preparing the X data ------------------

  # remove excessive variables
  X_imp <- cleanup(X_imp)

  # standardise the covariates in X (which are numeric and no intercept)
  X_imp_stand <- stand(X_imp)


  n <- length(y_imp)
  #define a place holder (ph)
  ph <- sample_imp(y_imp)[, 1]

  lmstart_all <- stats::lm(ph ~ 0 +., data = X_imp)

  X_model_matrix_1_all <- stats::model.matrix(lmstart_all)
  xnames_1 <- paste("X", 1:ncol(X_model_matrix_1_all), sep = "")

  unneeded <- is.na(stats::coefficients(lmstart_all))
  xnames_2 <- xnames_1[!unneeded]

  tmp_2_all <- data.frame(y = y_imp)
  tmp_2_all[, xnames_2] <- X_model_matrix_1_all[, !unneeded, drop = FALSE]

  tmp_2_sub <- tmp_2_all[!missind, , drop = FALSE]

  everything <- mice::mice(data = tmp_2_all, m = 1,
                     method = "norm",
                     predictorMatrix = (1 - diag(1, ncol(tmp_2_all))),
                     visitSequence = (1:ncol(tmp_2_all))[apply(is.na(tmp_2_all),2,any)],
                     post = vector("character", length = ncol(tmp_2_all)),
                     defaultMethod = "norm",
                     maxit = 10,
                     diagnostics = TRUE,
                     printFlag = FALSE,
                     seed = NA,
                     imputationMethod = NULL,
                     defaultImputationMethod = NULL,
                     data.init = NULL)


  y_ret <- data.frame(y_ret = y_imp)
  y_ret[missind, 1] <- everything$imp[[1]][, 1]

  return(y_ret)

}


# Generate documentation with devtools::document()
# Build package with devtools::build() and devtools::build(binary = TRUE) for zips
