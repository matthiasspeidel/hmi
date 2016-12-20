#################################################
#########Imputation function#####################
#################################################



#' The function for hierarchical imputation of contious variables.
#'
#' The function is called by the wrapper.
#' @param y_imp_multi A Vector with the variable to impute.
#' @param X_imp_multi A data.frame with the fixed effects variables.
#' @param M An integer defining the number of imputations that should be made.
#' @return A n x M matrix. Each column is one of M imputed y-variables.
imp_cont <- function(y_imp_multi,
                      X_imp_multi,
                      M = 10){

  # -----------------------------preparing the data ------------------
  # -- standardise the covariates in X (which are numeric and no intercept)
  types <- array(dim = ncol(X_imp_multi))
  for(j in 1:length(types)) types[j] <- get_type(X_imp_multi[, j])
  need_stand <- types == "cont"

  X_imp_multi_stand <- X_imp_multi
  X_imp_multi_stand[, need_stand] <- scale(X_imp_multi[, need_stand, drop = FALSE])
  #X_imp_multi %>% mutate_each_(funs(scale), vars = names(need_stand)[need_stand])
  #generate model.matrix (from the class matrix)
  n <- nrow(X_imp_multi_stand)
  tmp_0 <- data.frame(y_binary = rnorm(n), X_imp_multi_stand)
  X_model_matrix_1 <- model.matrix(y_binary ~ 0 +., data = tmp_0)

  # Remove ` from the variable names
  colnames(X_model_matrix_1) <- gsub("`", "", colnames(X_model_matrix_1))


  # -------------- calling the gibbs sampler to get imputation parameters----

  tmp_1 <- data.frame(y = y_imp_multi)
  xnames_1 <- paste("X", 1:ncol(X_model_matrix_1), sep = "")
  tmp_1[, xnames_1] <- X_model_matrix_1


  fixformula_1 <- formula(paste("y ~ 0 +", paste(xnames_1, collapse = "+"), sep = ""))

  reg_1 <- lm(fixformula_1, data = tmp_1, x = TRUE)
  # remove variables mit NA regression coefficient
  tmp_2 <- data.frame(y = y_imp_multi)
  xnames_2 <- xnames_1[!is.na(coefficients(reg_1))]
  tmp_2[, xnames_2] <- X_model_matrix_1[, !is.na(coefficients(reg_1)), drop = FALSE]


  fixformula_2 <- formula(paste("y ~ 0 +", paste(xnames_2, collapse = "+"), sep = ""))

  reg_2 <- lm(fixformula_2, data = tmp_2, x = TRUE)
  X_model_matrix_2 <- model.matrix(reg_2)


  reg_xtx <- solve(t(reg_2$x) %*% reg_2$x)

  dfs <- reg_2$df.residual

  reg_var <- sigma(reg_2)

  reg_coef <- coefficients(reg_2)
 # reg_coef_var <- vcov(reg)

  # -------------------- drawing samples with the parameters from the gibbs sampler --------
  y_imp <- array(NA, dim = c(n, M))
  ###start imputation
  for (j in 1:M){

    newsigma <- dfs * reg_var/rchisq(n = 1, dfs)

    fix_eff_imp <- MASS::mvrnorm(1, reg_coef, reg_xtx * newsigma)



    y_temp <- rnorm(n, X_model_matrix_2 %*% fix_eff_imp, sd = newsigma)

    y_imp[, j] <- ifelse(is.na(y_imp_multi), y_temp, y_imp_multi)
  }

  # --------- returning the imputed data --------------
  return(y_imp)
}


# Generate documentation with devtools::document()
# Build package with devtools::build() and devtools::build(binary = TRUE) for zips
