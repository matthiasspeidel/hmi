#' The function to impute ordered categorical variables
#'
#' The function uses the proportional odds logistic regression (polr) approach,
#' implemented in \code{mice}.
#' @param y_imp A Vector with the variable to impute.
#' @param X_imp A data.frame with the fixed effects variables.
#' @return A n x 1 data.frame with the original and imputed values as a factor.
imp_orderedcat_single <- function(y_imp, X_imp){


  categories <- levels(y_imp)


  # ----------------------------- preparing the X data ------------------
  # remove excessive variables
  X_imp <- cleanup(X_imp)

  # standardize X
  X_imp_stand <- stand(X_imp)

  #the missing indactor indicates, which values of y are missing.
  missind <- is.na(y_imp)

  #starting model
  ph <- sample_imp(y_imp)[, 1]
  tmp_0_all <- data.frame(target = ph)
  xnames_0 <- paste("X", 1:ncol(X_imp_stand), sep = "")
  tmp_0_all[xnames_0] <- X_imp_stand
  tmp_0_sub <- tmp_0_all[!missind, , drop = FALSE]

  reg_1_all <- MASS::polr(target ~ 1 + ., data =  tmp_0_all, method = "probit")
  reg_1_sub <- MASS::polr(target ~ 1 + ., data =  tmp_0_sub, method = "probit")

  X_model_matrix_1_all <- stats::model.matrix(reg_1_all)
  xnames_1 <- paste("X", 1:ncol(X_model_matrix_1_all), sep = "")


  #remove unneeded variables
  xnames_2 <- xnames_1[!is.na(stats::coefficients(reg_1_sub))]

  tmp_2_all <- data.frame(target = as.factor(y_imp)) # mice needs the variable as a factor
  tmp_2_all[, xnames_2] <- X_model_matrix_1_all[, !is.na(stats::coefficients(reg_1_sub)),
                                                drop = FALSE]

  everything <- mice::mice(data = tmp_2_all, m = 1,
              method = "polr",
              predictorMatrix = (1 - diag(1, ncol(tmp_2_all))),
              visitSequence = (1:ncol(tmp_2_all))[apply(is.na(tmp_2_all),2,any)],
              post = vector("character", length = ncol(tmp_2_all)),
              defaultMethod = "polr",
              maxit = 10,
              diagnostics = TRUE,
              printFlag = FALSE,
              seed = NA,
              imputationMethod = NULL,
              defaultImputationMethod = NULL,
              data.init = NULL)


  #Initialising the returning vector
  y_ret <- data.frame(y_ret = y_imp)

  y_ret[missind, 1] <- everything$imp[[1]][, 1]

  return(y_ret)
}




