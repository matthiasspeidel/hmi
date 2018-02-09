#' The function to impute ordered categorical variables
#'
#' The function uses the proportional odds logistic regression (polr) approach,
#' implemented in \code{mice}.
#' @param y_imp A Vector with the variable to impute.
#' @param X_imp A data.frame with the fixed effects variables.
#' @param pvalue A numeric between 0 and 1 denoting the threshold of p-values a variable in the imputation
#' model should not exceed. If they do, they are excluded from the imputation model.
#' @param rounding_degrees A numeric vector with the presumed rounding degrees.
#' @return A n x 1 data.frame with the original and imputed values as a factor.
imp_orderedcat_single <- function(y_imp,
                                  X_imp,
                                  pvalue = 0.2,
                                  rounding_degrees = c(1, 10, 100, 1000)){


  categories <- levels(y_imp)


  # ----------------------------- preparing the X data ------------------
  # remove excessive variables
  X_imp <- cleanup(X_imp)

  # standardize X
  X_imp_stand <- stand(X_imp, rounding_degrees = rounding_degrees)

  intercept_indicator <- apply(X_imp_stand, 2, get_type) == "intercept"
  global_intercept_indicator <- intercept_indicator
  X_imp_stand <- X_imp_stand[, !intercept_indicator, drop = FALSE]
  #the missing indactor indicates, which values of y are missing.
  missind <- is.na(y_imp)

  n <- length(y_imp)

  # ----------- set up a maximal model matrix with all possible relevant (dummy) variables -----
  # In the imputation model only actually relevant (dummy) variables shall be present.
  # THis is done by setting up a mirror of the initial model matrix.
  # Then step by step this model matrix is reduced to all actually relevant (dummy) variables.
  # This reduction is based on models using the observed data.
  # The last step prior to the imputation-parameters estimation is to restrict the initial mode matrix
  # to those variables, left in the reduced mirror model matrix.
  #define a place holder (ph)
  ph <- sample_imp(y_imp)[, 1]

  tmp_0_all <- data.frame(target = ph, X_imp_stand)
  xnames_1 <- colnames(X_imp_stand)

  if(any(intercept_indicator)){
    tmp_formula <- paste("target~ 1 + ", paste(xnames_1, collapse = "+"), sep = "")
  }else{
    tmp_formula <- paste("target~ 0 + ", paste(xnames_1, collapse = "+"), sep = "")
  }

  oldw <- getOption("warn")
  options(warn = -1)
  reg_1_all <- MASS::polr(stats::formula(tmp_formula), data = tmp_0_all, method = "probit")
  options(warn = oldw)

  X_model_matrix_1_all <- stats::model.matrix(reg_1_all)
  xnames_1 <- paste("X", 1:ncol(X_model_matrix_1_all), sep = "")
  colnames(X_model_matrix_1_all) <- xnames_1

  tmp_0_all <- data.frame(target = ph)
  tmp_0_all[, xnames_1] <- X_model_matrix_1_all

  #From this initial model matrix X_model_matrix_1_all
  #now step by step irrelavant variables are removed.
  X_model_matrix_1_sub <- X_model_matrix_1_all[!missind, , drop = FALSE]

  # The first step of the reduction is to remove variables having a non-measurable effect
  # (e.g. due to colinearity) on y.
  # tmp_1 shall include the covariates (like X_model_matrix) and additionally the target variable
  ph_sub <- ph[!missind]
  tmp_1_sub <- data.frame(target = ph_sub)

  intercept_indicator <- apply(X_model_matrix_1_sub, 2, get_type) == "intercept"
  X_model_matrix_1_sub <- X_model_matrix_1_sub[, !intercept_indicator, drop = FALSE]
  xnames_1 <- colnames(X_model_matrix_1_sub)
  tmp_1_sub[, xnames_1] <- X_model_matrix_1_sub

  if(any(intercept_indicator)){
    tmp_formula <- paste("target~ 1 + ", paste(xnames_1, collapse = "+"), sep = "")
  }else{
    tmp_formula <- paste("target~ 0 + ", paste(xnames_1, collapse = "+"), sep = "")
  }

  oldw <- getOption("warn")
  options(warn = -1)
  reg_1_sub <- MASS::polr(stats::formula(tmp_formula), data = tmp_1_sub, method = "probit")
  options(warn = oldw)
  #remove unneeded variables
  X_model_matrix_1_sub <- X_model_matrix_1_sub[, !is.na(stats::coefficients(reg_1_sub)),
                                               drop = FALSE]
  ############################################################
  # Remove insignificant variables from the imputation model #
  ############################################################
  check <- TRUE
  while(check){
    tmp_1_sub <- data.frame(target = ph_sub)

    X_model_matrix_1_sub <- X_model_matrix_1_sub[,
                                    apply(X_model_matrix_1_sub, 2, get_type) != "intercept",
                                    drop = FALSE]
    xnames_1 <- colnames(X_model_matrix_1_sub)
    tmp_1_sub[, xnames_1] <- X_model_matrix_1_sub

    if(any(global_intercept_indicator)){
      tmp_formula <- paste("target ~ 1 + ", paste(xnames_1, collapse = "+"), sep = "")
    }else{
      tmp_formula <- paste("target ~ 0 + ", paste(xnames_1, collapse = "+"), sep = "")
    }

    oldw <- getOption("warn")
    options(warn = -1)
    reg_1_sub <- MASS::polr(stats::formula(tmp_formula), data = tmp_1_sub, method = "probit",
                            model = TRUE, Hess = TRUE)
    options(warn = oldw)

    #stats::pnorm(abs(stats::coef(summary(reg_1_sub)))[, "t value"], lower.tail = FALSE) * 2
    tmp <- abs(stats::coef(summary(reg_1_sub))[, "t value"])
    #Evaluate only the real covariates that went into the model, but not the newly generated variables
    #like "A|B" and "B|C" (in the case the target variable consited of the ordered categories "A"<"B"<"C")
    pvalues <- stats::pt(tmp[names(tmp) %in% xnames_1], df = stats::df.residual(reg_1_sub),
                         lower.tail = FALSE) * 2
    insignificant_variables <- which(pvalues > pvalue)
    most_insignificant <- insignificant_variables[which.max(pvalues[insignificant_variables])]

    if(length(most_insignificant) == 0){
      check <- FALSE
    }else{
      tmp <- stats::model.matrix(reg_1_sub) #if an additional intercept variable is included by the model
      #we cannot run stats::model.matrix(reg_1_sub)[, -most_insignificant]
      #Because most_insignificant refers to a situation without an intercept variable.
      X_model_matrix_1_sub <- tmp[, !colnames(tmp) %in% names(most_insignificant), drop = FALSE]
    }
  }

  tmp_2_all <- tmp_0_all[, colnames(tmp_1_sub), drop = FALSE]
  tmp_2_all$target[missind] <- NA

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




