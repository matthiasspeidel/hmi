#' The function for imputation of binary variables.
#'
#' The function is called by the wrapper.
#' @param y_imp A Vector with the variable to impute.
#' @param X_imp A data.frame with the fixed effects variables.
#' @return A n x 1 data.frame with the original and imputed values.
imp_binary_single <- function(y_imp,
                      X_imp){


  # ----------------------------- preparing the y data ------------------
  # stransform y_imp into a real binary with only zeros and ones (and NAs).
  first_possibility <- utils::head(sort(y_imp), n = 1)
  second_possibility <- utils::tail(sort(y_imp), n = 1)
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
  X_imp <- remove_excessives(X_imp)

  # standardize X
  X_imp_stand <- stand(X_imp)

  #the missing indactor indicates, which values of y are missing.
  missind <- is.na(y_binary)
  n <- length(y_binary)

  #define a place holder (ph)
  ph <- sample_imp(y_binary[, 1])[, 1]

  tmp_0_sub <- data.frame(target = y_binary, X_imp_stand)[!missind, , drop = FALSE]
  tmp_0_all <- data.frame(target = ph, X_imp_stand)

  X_model_matrix_1_sub <- stats::model.matrix(target ~ 0 + ., data = tmp_0_sub)
  X_model_matrix_1_all <- stats::model.matrix(target ~ 0 + ., data = tmp_0_all)

  # Remove ` from the variable names
  colnames(X_model_matrix_1_sub) <- gsub("`", "", colnames(X_model_matrix_1_sub))
  colnames(X_model_matrix_1_all) <- gsub("`", "", colnames(X_model_matrix_1_all))

  xnames_1 <- paste("X", 1:ncol(X_model_matrix_1_sub), sep = "")

  #tmp_1 shall include the covariates (like X_model_matrix) and additionally the target variable
  tmp_1_all <- data.frame(target = y_binary)
  tmp_1_all[, xnames_1] <- X_model_matrix_1_all
  tmp_1_sub <- tmp_1_all[!missind, , drop = FALSE]


  glm_fixpart_1 <- paste("target~ 0 + ", paste(xnames_1, collapse = "+"), sep = "")

  reg_1_sub <- stats::glm(stats::formula(glm_fixpart_1), data = tmp_1_sub,
                          family = stats::binomial(link = "logit"))

  #remove unneeded variables

  xnames_2 <- xnames_1[!is.na(stats::coefficients(reg_1_sub))]

  tmp_2_all <- data.frame(target = y_binary)
  tmp_2_all[, xnames_2] <- X_model_matrix_1_all[, !is.na(stats::coefficients(reg_1_sub)),
                                                drop = FALSE]

  everything <- mice::mice(data = tmp_2_all, m = 1,
                     method = "logreg",
                     predictorMatrix = (1 - diag(1, ncol(tmp_2_all))),
                     visitSequence = (1:ncol(tmp_2_all))[apply(is.na(tmp_2_all),2,any)],
                     post = vector("character", length = ncol(tmp_2_all)),
                     defaultMethod = "logreg",
                     maxit = 10,
                     diagnostics = TRUE,
                     printFlag = FALSE,
                     seed = NA,
                     imputationMethod = NULL,
                     defaultImputationMethod = NULL,
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

  return(y_ret)
}


# Generate documentation with devtools::document()
# Build package with devtools::build() and devtools::build(binary = TRUE) for zips
