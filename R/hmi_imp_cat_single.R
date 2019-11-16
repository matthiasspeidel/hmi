#' The function to impute unordered categorical variables
#'
#' The function uses regression trees for imputation implemented in \code{mice}.
#' The principle is the following:
#' For each observation it is calculated at which leave it would end.
#' Then one (randomly selected) observation of the other observations found on this leave
#' functions as a donor.
#' @param y_imp A Vector with the variable to impute.
#' @param X_imp A data.frame with the fixed effects variables.
#' @param pvalue A numeric between 0 and 1 denoting the threshold of p-values a variable in the imputation
#' model should not exceed. If they do, they are excluded from the imputation model.
#' @param k An integer defining the allowed maximum of levels in a factor covariate.
#' @return A n x 1 data.frame with the original and imputed values.
imp_cat_single <- function(y_imp,
                           X_imp,
                           pvalue = 0.2,
                           k = Inf){

  if(min(table(y_imp)) < 2) {
    stop("Too few observations per category in a categorical target variable.")
  }

  # ----------------------------- preparing the X data ------------------
  # remove excessive variables
  X_imp <- cleanup(X_imp, k = k)

  # standardize X
  X_imp_stand <- stand(X_imp)


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

  tmp_formula <- paste("target~ 0 + ", paste(xnames_1, collapse = "+"), sep = "")
  reg_1_all <- nnet::multinom(stats::formula(tmp_formula), data = tmp_0_all, trace = FALSE)

  X_model_matrix_1_all <- stats::model.matrix(reg_1_all)
  xnames_1 <- paste("X", 1:ncol(X_model_matrix_1_all), sep = "")
  colnames(X_model_matrix_1_all) <- xnames_1

  tmp_0_all <- data.frame(target = ph)
  tmp_0_all[, xnames_1] <- as.data.frame(X_model_matrix_1_all)

  #From this initial model matrix X_model_matrix_1_all
  #now step by step irrelavant variables are removed.
  X_model_matrix_1_sub <- X_model_matrix_1_all[!missind, , drop = FALSE]

  # The first step of the reduction is to remove variables having a non-measurable effect
  # (e.g. due to colinearity) on y.
  # tmp_1 shall include the covariates (like X_model_matrix) and additionally the target variable
  ph_sub <- ph[!missind]
  tmp_1_sub <- data.frame(target = ph_sub)
  xnames_1 <- colnames(X_model_matrix_1_sub)
  tmp_1_sub[, xnames_1] <- as.data.frame(X_model_matrix_1_sub)

  tmp_formula <- paste("target~ 0 + ", paste(xnames_1, collapse = "+"), sep = "")
  reg_1_sub <- nnet::multinom(stats::formula(tmp_formula), data = tmp_1_sub, trace = FALSE)

  #remove unneeded variables
  tmp <- stats::coefficients(reg_1_sub)
  X_model_matrix_1_sub <- X_model_matrix_1_sub[, !apply(tmp, 2, function(x) any(is.na(x))),
                                               drop = FALSE]

  # Remove insignificant variables from the imputation model
  check <- TRUE
  while(check){
    tmp_1_sub <- data.frame(target = ph_sub)
    xnames_1 <- colnames(X_model_matrix_1_sub)
    tmp_1_sub[, xnames_1] <- X_model_matrix_1_sub
    tmp_formula <- paste("target~ 0 + ", paste(xnames_1, collapse = "+"), sep = "")
    reg_1_sub <- nnet::multinom(stats::formula(tmp_formula), data = tmp_1_sub, trace = FALSE)

    z <- summary(reg_1_sub)$coefficients / summary(reg_1_sub)$standard.errors
    pvalues <- apply((1 - stats::pnorm(abs(z)))*2, 2, min)
    insignificant_variables <- which(pvalues > pvalue)
    most_insignificant <- insignificant_variables[which.max(pvalues[insignificant_variables])]

    if(length(most_insignificant) == 0){
      check <- FALSE
    }else{
      tmp <- stats::model.matrix(reg_1_sub) #if an additional intercept variable is included by the model
      #we cannot run stats::model.matrix(reg_1_sub)[, -most_insignificant]
      #Because most_insignificant refers to a situation without an intercept variable.

      #.. drop the insignificant variable from the model.matrix, but only if at least 1 variable remains
      tmp_MM <- tmp[, !colnames(tmp) %in% names(most_insignificant), drop = FALSE]
      if(ncol(tmp_MM) == 0){
        check <- FALSE
      }else{
        X_model_matrix_1_sub <- tmp_MM
      }
    }
  }

  tmp_2_all <- tmp_0_all[, colnames(tmp_1_sub), drop = FALSE]
  tmp_2_all$target[missind] <- NA

  #Constant variables are of no use in cart models. Therefore they are removed.
  tmp_2_all <- tmp_2_all[, get_type(tmp_2_all) != "intercept", drop = FALSE]

  #Initialising the returning vector
  y_ret <- y_imp

  #if no covariates are left, a sample imputation is needed.
  if(ncol(tmp_2_all) < 2){
    y_ret <- sample_imp(y_imp)[, 1]
  }else{

    # run cart imputation from mice.
    # Note: mice needs the variable to impute as a factor.
    # Therefore if the data passed to hmi weren't factors (but a vector characters),
    # the factors returned by mice need to be translated back into characters.
    # In case a factor was passed to hmi, it still works to pass a character to this factor. cf
    # a <- factor(c("x", "y", "z"))
    # a[2] <- "z"
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
                             data.init = NULL)
    y_ret[missind] <- as.character(everything$imp[[1]][, 1])
  }



  ret <- data.frame(y_ret = y_ret, stringsAsFactors = FALSE)
  return(ret)
}
