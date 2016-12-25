#' The function to impute categorical variables
#'
#' The function uses regression trees for imputation. The principle is the following:
#' For each observation it is calculated at wich leave it would end.
#' Then one (randomly selected) observation of the other observations found on this leave,
#' functions as a donor.
#' @param y_imp_multi A Vector with the variable to impute.
#' @param X_imp_multi A data.frame with the fixed effects variables.
#' @return A n x 1 data.frame
imp_cat_cart <- function(y_imp_multi,
                         X_imp_multi){
  #Initialising the returning vector
  y_imp <- as.matrix(y_imp_multi, ncol = 1)

  #the missing indactor indicates, which values of y are missing.
  missind <- is.na(y_imp_multi)


  types <- array(dim = ncol(X_imp_multi))
  for(j in 1:length(types)) types[j] <- get_type(X_imp_multi[, j])
  #need_stand <- types == "cont"
  categorical <- types == "categorical"

  #remove categories with more than 10 observations as the model in the current form
  #will cause later numerical probles
  too_many_levels <- colnames(X_imp_multi[, categorical, drop = FALSE])[
    apply(X_imp_multi[, categorical, drop = FALSE], 2, function(x) nlevels(factor(x))) > 10]
  X_imp_multi <- X_imp_multi[, !names(X_imp_multi) %in% too_many_levels, drop = FALSE]


  n <- length(y_imp_multi)
  lmstart <- lm(rnorm(n) ~ 0 +., data = X_imp_multi)

  X_model_matrix_1 <- model.matrix(lmstart)
  xnames_1 <- paste("X", 1:ncol(X_model_matrix_1), sep = "")

  tmp_1 <- data.frame(y = rnorm(n))
  tmp_1[, xnames_1] <- X_model_matrix_1

  reg_1 <- lm(y ~ 0 + . , data = tmp_1)


  tmp_2 <- data.frame(y = as.factor(y_imp_multi)) # mice needs the binary variable as a factor

  xnames_2 <- xnames_1[!is.na(coefficients(reg_1))]
  tmp_2[, xnames_2] <- X_model_matrix_1[, !is.na(coefficients(reg_1)), drop = FALSE]

  everything <- mice(data = tmp_2, m = 1,
                     method = "cart",
                     predictorMatrix = (1 - diag(1, ncol(tmp_2))),
                     visitSequence = (1:ncol(tmp_2))[apply(is.na(tmp_2),2,any)],
                     post = vector("character", length = ncol(tmp_2)),
                     defaultMethod = "cart",
                     maxit = 10,
                     diagnostics = TRUE,
                     printFlag = FALSE,
                     seed = NA,
                     imputationMethod = NULL,
                     defaultImputationMethod = NULL,
                     data.init = NULL)

  y_imp[is.na(y_imp_multi), 1] <- everything$imp[[1]][, 1]

  return(data.frame(y_imp = factor(y_imp)))
}
