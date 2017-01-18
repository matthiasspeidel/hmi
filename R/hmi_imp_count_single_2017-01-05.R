#' The function for imputation of binary variables.
#'
#' The function is called by the wrapper.
#' @param y_imp_multi A Vector with the variable to impute.
#' @param X_imp_multi A data.frame with the fixed effects variables.
#' @return A n x 1 matrix. The column is one set of imputed y-variables.
imp_count_single <- function(y_imp_multi,
                       X_imp_multi){

  need_stand <- apply(X_imp_multi, 2, get_type) == "cont"
  X_imp_multi_stand <- X_imp_multi
  X_imp_multi_stand[, need_stand] <- scale(X_imp_multi[, need_stand])

  #generate model.matrix (from the class matrix)
  n <- nrow(X_imp_multi_stand)
  X_model_matrix <- stats::model.matrix(stats::rnorm(n) ~ 0 + ., data = X_imp_multi_stand)
  # Remove ` from the variable names
  colnames(X_model_matrix) <- gsub("`", "", colnames(X_model_matrix))

  # -------------- calling the gibbs sampler to get imputation parameters----


  #n <- length(y_imp_multi)
  lmstart <- stats::lm(stats::rnorm(n) ~ 0 +., data = X_imp_multi)

  X_model_matrix_1 <- stats::model.matrix(lmstart)
  xnames_1 <- paste("X", 1:ncol(X_model_matrix_1), sep = "")

  tmp_1 <- data.frame(y = stats::rnorm(n))
  tmp_1[, xnames_1] <- X_model_matrix_1

  reg_1 <- stats::lm(y ~ 0 + . , data = tmp_1)

  blob <- y_imp_multi
  tmp_2 <- data.frame(target = blob)

  xnames_2 <- xnames_1[!is.na(stats::coefficients(reg_1))]
  X_model_matrix_2 <- X_model_matrix_1[, !is.na(stats::coefficients(reg_1)), drop = FALSE]
  tmp_2[, xnames_2] <- X_model_matrix_2

  fixformula <- stats::formula(paste("target~", paste(xnames_2, collapse = "+"), "- 1", sep = ""))

  prior <- list(R = list(V = 1e-07, nu = -2))

  MCMCglmm_draws <- MCMCglmm::MCMCglmm(fixformula, data = tmp_2,
                                       verbose = FALSE, pr = TRUE, prior = prior,
                                       family = "poisson",
                                       saveX = TRUE,
                                       nitt = 3000,
                                       thin = 10,
                                       burnin = 1000)

  pointdraws <- MCMCglmm_draws$Sol
  xdraws <- pointdraws[, 1:ncol(X_model_matrix_2), drop = FALSE]
  variancedraws <- MCMCglmm_draws$VCV
  # the last column contains the variance (not standard deviation) of the residuals

  number_of_draws <- nrow(pointdraws)
  select.record <- sample(1:number_of_draws, 1, replace = TRUE)

  # -------------------- drawing samples with the parameters from the gibbs sampler --------
  ###start imputation


  fix.eff.imp <- matrix(xdraws[select.record, ], nrow = ncol(X_model_matrix_2))

  sigma.y.imp <- sqrt(variancedraws[select.record, ncol(variancedraws)])


  lambda <- exp(stats::rnorm(n, X_model_matrix_2 %*% fix.eff.imp, sigma.y.imp))

  # note: the maximum lambda is something about 2.14e9
  if(max(lambda) > 2.14e9) {
    stop("Imputation of count variable failed due to a too high value of lambda.
This can occur if an observation in your data is an outlier regarding the covariates of the imputation model.
What again can be caused by highly variant imputation models for these covariates due high missing rates.")
  }
  y_imp <- ifelse(is.na(y_imp_multi), stats::rpois(n, lambda), y_imp_multi)
  y_imp <- as.matrix(y_imp, ncol = 1)


  # --------- returning the imputed data --------------
  return(y_imp)
}



# Generate documentation with devtools::document()
# Build package with devtools::build() and devtools::build(binary = TRUE) for zips
