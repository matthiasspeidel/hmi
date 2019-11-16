#' The function for imputation of binary variables.
#'
#' The function is called by the wrapper.
#' @param y_imp A Vector with the variable to impute.
#' @param X_imp A data.frame with the fixed effects variables.
#' @param pvalue A numeric between 0 and 1 denoting the threshold of p-values a variable in the imputation
#' model should not exceed. If they do, they are excluded from the imputation model.
#' @param k An integer defining the allowed maximum of levels in a factor covariate.
#' @return A n x 1 data.frame with the original and imputed values.
imp_binary_single <- function(y_imp,
                      X_imp,
                      pvalue = 0.2,
                      k = Inf){


  # ----------------------------- preparing the y data ------------------
  # stransform y_imp into a real binary with only zeros and ones (and NAs).
  first_possibility <- sort(y_imp)[1]
  second_possibility <- rev(sort(y_imp))[1]
  y_binary <- data.frame(target = factor(y_imp, labels = c(0, 1)))

  # If one category has less then two observations, no binary model can be estimated.
  # So the imputation routines has to stop.
  if(min(table(y_binary)) < 2){
    stop("A binary (or maybe a semicontinuous) variable has less than two observations in one category.
         Consider removing this variable
         (or in the case of a semicontinuous variable, to specify it as continuous in the list_of_types (see ?hmi)).")
  }

  # ----------------------------- preparing the X data ------------------
  # remove excessive variables
  X_imp <- cleanup(X_imp, k = k)

  # standardize X
  X <- stand(X_imp)

  #the missing indactor indicates, which values of y are missing.
  missind <- is.na(y_binary)
  n <- length(y_binary)

  # ----------- set up a maximal model matrix with all possible relevant (dummy) variables -----
  # In the imputation model only actually relevant (dummy) variables shall be present.
  # THis is done by setting up a mirror of the initial model matrix.
  # Then step by step this model matrix is reduced to all actually relevant (dummy) variables.
  # This reduction is based on models using the observed data.
  # The last step prior to the imputation-parameters estimation is to restrict the initial mode matrix
  # to those variables, left in the reduced mirror model matrix.

  #define a place holder (ph)
  ph <- sample_imp(y_binary[, 1])[, 1]

  tmp_0_all <- data.frame(target = ph, X)
  xnames_1 <- colnames(X)

  tmp_formula <- paste("target~ 0 + ", paste(xnames_1, collapse = "+"), sep = "")
  reg_1_all <- stats::glm(stats::formula(tmp_formula), data = tmp_0_all,
                          family = stats::binomial(link = "logit"))

  X_model_matrix_1_all <- stats::model.matrix(reg_1_all)
  xnames_1 <- paste("X", 1:ncol(X_model_matrix_1_all), sep = "")
  colnames(X_model_matrix_1_all) <- xnames_1

  tmp_0_all <- data.frame(target = ph)
  tmp_0_all[, xnames_1] <- as.data.frame(X_model_matrix_1_all)

  #From this initial model matrix X_model_matrix_1_all
  #now step by step irrelavant variables are removed.
  X_model_matrix_1_sub <- X_model_matrix_1_all[!missind, , drop = FALSE]


  #first step of the reduction is to remove variables having a NA-effect (e.g. due to colinearity) on y
  #tmp_1 shall include the covariates (like X_model_matrix) and additionally the target variable
  ph_sub <- y_binary[!missind, , drop = FALSE]
  tmp_1_sub <- data.frame(target = ph_sub)
  xnames_1 <- colnames(X_model_matrix_1_sub)
  tmp_1_sub[, xnames_1] <- as.data.frame(X_model_matrix_1_sub)

  tmp_formula <- paste("target~ 0 + ", paste(xnames_1, collapse = "+"), sep = "")

  reg_1_sub <- stats::glm(stats::formula(tmp_formula), data = tmp_1_sub,
                          family = stats::binomial(link = "logit"))

  #remove unneeded variables
  X_model_matrix_1_sub <- X_model_matrix_1_sub[, !is.na(stats::coefficients(reg_1_sub)),
                       drop = FALSE]

  # Remove insignificant variables from the imputation model
  check <- TRUE
  while(check){
    tmp_1_sub <- data.frame(target = ph_sub)
    xnames_1 <- colnames(X_model_matrix_1_sub)
    tmp_1_sub[, xnames_1] <- X_model_matrix_1_sub
    tmp_formula <- paste("target~ 0 + ", paste(xnames_1, collapse = "+"), sep = "")
    reg_1_sub <- stats::glm(stats::formula(tmp_formula), data = tmp_1_sub,
                              family = stats::binomial(link = "logit"))

    pvalues <- summary(reg_1_sub)$coefficients[, 4]
    insignificant_variables <- which(pvalues > pvalue)
    most_insignificant <- insignificant_variables[which.max(pvalues[insignificant_variables])]

    if(length(most_insignificant) == 0){
      check <- FALSE
    }else{
      #.. drop the insignificant variable from the model.matrix, but only if at least 1 variable remains
      tmp_MM <- stats::model.matrix(reg_1_sub)[, -most_insignificant, drop = FALSE]
      if(ncol(tmp_MM) == 0){
        check <- FALSE
      }else{
        X_model_matrix_1_sub <- tmp_MM
      }
    }

  }

  tmp_2_all <- tmp_0_all[, colnames(tmp_1_sub), drop = FALSE]
  tmp_2_all$target[missind] <- NA

  everything <- mice::mice(data = tmp_2_all,
                           m = 1,
                           method = "logreg",
                           predictorMatrix = (1 - diag(1, ncol(tmp_2_all))),
                           visitSequence = (1:ncol(tmp_2_all))[apply(is.na(tmp_2_all), 2, any)],
                           post = vector("character", length = ncol(tmp_2_all)),
                           defaultMethod = "logreg",
                           maxit = 10,
                           diagnostics = TRUE,
                           printFlag = FALSE,
                           seed = NA,
                           data.init = NULL)

  indicator <- as.numeric(as.character(mice::complete(everything, 1)$target))

  #Initialising the returning vector
  y_ret <- as.data.frame(y_imp)

  #ifelse(indicator == 1, first_possibility, second_possibility) doesn't work with factors
  # in a way I would need it.
  for(i in which(is.na(y_imp))){
    if(indicator[i] == 0){
      y_ret[i, 1] <- first_possibility
    }else{
      y_ret[i, 1] <- second_possibility
    }
  }

  colnames(y_ret) <- "y_ret"
  return(y_ret)
}


# Generate documentation with devtools::document()
# Build package with devtools::build() and devtools::build(binary = TRUE) for zips
