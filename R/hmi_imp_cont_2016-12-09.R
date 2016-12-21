#################################################
#########Imputation function#####################
#################################################



#' The function for hierarchical imputation of contious variables.
#'
#' The function is called by the wrapper.
#' @param y_imp_multi A Vector with the variable to impute.
#' @param X_imp_multi A data.frame with the fixed effects variables.
#' @param M An integer defining the number of imputations that should be made.
#' @param allowed_max_value A single numeric Value which shall not be exceeded
#' when values are imputed (e.g. the age of a person can be limited to 125).
#' @param allowed_max_variable A character naming a variable V.
#' For each Y_i the value of V_i shall not exceeded
#' (e.g. the net income shall not exceed the gross income).
#' Note that a new imputed value has to satisfy both conditions of \code{allowed_max_value}
#' and \code{allowed_max_variable} at the same time.
#' @param allowed_min_value Analog to \code{allowed_max_value}.
#' @param allowed_min_variable Analog to \code{allowed_max_variable}.
#' @return A n x M matrix. Each column is one of M imputed y-variables.
imp_cont <- function(y_imp_multi,
                      X_imp_multi,
                      M = 10,
                      allowed_max_value = Inf,
                      allowed_max_variable = NULL,
                      allowed_min_value = -Inf,
                      allowed_min_variable = NULL){

  # -----------------------------preparing the data ------------------
  # -- standardise the covariates in X (which are numeric and no intercept)
  types <- array(dim = ncol(X_imp_multi))
  for(j in 1:length(types)) types[j] <- get_type(X_imp_multi[, j])
  need_stand <- types == "cont"

  X_imp_multi_stand <- X_imp_multi
  X_imp_multi_stand[, need_stand] <- scale(X_imp_multi[, need_stand])
  #X_imp_multi %>% mutate_each_(funs(scale), vars = names(need_stand)[need_stand])
  #generate model.matrix (from the class matrix)
  n <- nrow(X_imp_multi_stand)
  tmp0 <- data.frame(y_binary = rnorm(n), X_imp_multi_stand)
  X_model_matrix <- model.matrix(y_binary ~ 0 +., data = tmp0)
  #!!! BESSER LOESEN!!!
  #alt: model.matrix(rnorm(n) ~ 0 + ., data = X_imp_multi_stand)
  # Remove ` from the variable names
  colnames(X_model_matrix) <- gsub("`", "", colnames(X_model_matrix))


  # -------------- calling the gibbs sampler to get imputation parameters----

  tmp <- data.frame(y = y_imp_multi)
  xnames <- paste("X", 1:ncol(X_model_matrix), sep = "")
  tmp[, xnames] <- X_model_matrix


  fixformula <- formula(paste("y ~ 0 +", paste(xnames, collapse = "+"), sep = ""))

  reg <- lm(fixformula, data = tmp, x = TRUE)

  reg_xtx <- solve(t(reg$x) %*% reg$x)

  dfs <- reg$df.residual

  reg_var <- sigma(reg)

  reg_coef <- coefficients(reg)
 # reg_coef_var <- vcov(reg)

  # -------------------- drawing samples with the parameters from the gibbs sampler --------
  y_imp <- array(NA, dim = c(n, M))
  ###start imputation
  for (j in 1:M){

    newsigma <- dfs * reg_var/rchisq(n = 1, dfs)

    fix_eff_imp <- MASS::mvrnorm(1, reg_coef, reg_xtx * newsigma)



    y_temp <- rnorm(n, X_model_matrix %*% fix_eff_imp, sd = newsigma)

    y_imp[, j] <- ifelse(is.na(y_imp_multi), y_temp, y_imp_multi)
  }

  # --------- returning the imputed data --------------
  return(y_imp)
}


# Generate documentation with devtools::document()
# Build package with devtools::build() and devtools::build(binary = TRUE) for zips
