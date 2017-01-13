#' The function for imputation of binary variables.
#'
#' The function is called by the wrapper.
#' @param y_imp_multi A Vector with the variable to impute.
#' @param X_imp_multi A data.frame with the fixed effects variables.
#' @return A n x 1 data.frame The column is one set of imputed y-variables.
imp_binary_single <- function(y_imp_multi,
                      X_imp_multi){

  #Initialising the returning vector
  y_imp <- data.frame(y_imp_multi)

  #the missing indactor indicates, which values of y are missing.
  missind <- is.na(y_imp_multi)


  types <- array(dim = ncol(X_imp_multi))
  for(j in 1:length(types)) types[j] <- get_type(X_imp_multi[, j])

  categorical <- types == "categorical"

  #remove categories with more than 10 observations as the model in the current form
  #will cause later numerical probles
  too_many_levels <- colnames(X_imp_multi[, categorical, drop = FALSE])[
    apply(X_imp_multi[, categorical, drop = FALSE], 2, function(x) nlevels(factor(x))) > 10]
  X_imp_multi <- X_imp_multi[, !names(X_imp_multi) %in% too_many_levels, drop = FALSE]


  n <- length(y_imp_multi)
  lmstart <- stats::lm(stats::rnorm(n) ~ 0 + ., data = X_imp_multi)

  X_model_matrix_1 <- stats::model.matrix(lmstart)
  xnames_1 <- paste("X", 1:ncol(X_model_matrix_1), sep = "")

  tmp_1 <- data.frame(y = stats::rnorm(n))
  tmp_1[, xnames_1] <- X_model_matrix_1

  reg_1 <- stats::lm(y ~ 0 + . , data = tmp_1)

  # mice needs the binary variable as a factor
  # so we force the date to be (at least temporarily) factors
  # afterwards, we return the data into their original format
  first_possibility <- utils::head(sort(y_imp_multi), n = 1)
  second_possibility <- utils::tail(sort(y_imp_multi), n = 1)
  tmp_2 <- data.frame(y = factor(y_imp_multi, labels = c(1, 2)))



  xnames_2 <- xnames_1[!is.na(stats::coefficients(reg_1))]
  tmp_2[, xnames_2] <- X_model_matrix_1[, !is.na(stats::coefficients(reg_1)), drop = FALSE]

  everything <- mice::mice(data = tmp_2, m = 1,
                     method = "logreg",
                     predictorMatrix = (1 - diag(1, ncol(tmp_2))),
                     visitSequence = (1:ncol(tmp_2))[apply(is.na(tmp_2),2,any)],
                     post = vector("character", length = ncol(tmp_2)),
                     defaultMethod = "logreg",
                     maxit = 10,
                     diagnostics = TRUE,
                     printFlag = FALSE,
                     seed = NA,
                     imputationMethod = NULL,
                     defaultImputationMethod = NULL,
                     data.init = NULL)

  indicator <- as.numeric(as.character(mice::complete(everything, 1)$y))

  #ifelse(indicator == 1, first_possibility, second_possibility) doesn't work with factors
  # in a way I would need it.
  for(i in which(is.na(y_imp_multi))){
    if(indicator[i] == 1){
      y_imp[i, 1] <- first_possibility
    }else{
      y_imp[i, 1] <- second_possibility
    }
  }

  return(y_imp)
}


# Generate documentation with devtools::document()
# Build package with devtools::build() and devtools::build(binary = TRUE) for zips
