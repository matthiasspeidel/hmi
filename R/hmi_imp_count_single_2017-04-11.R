#' The function for imputation of binary variables.
#'
#' The function is called by the wrapper.
#' @param y_imp A vector with the variable to impute.
#' @param X_imp A data.frame with the fixed effects variables.
#' @return A n x 1 data.frame with the original and imputed values.
imp_count_single <- function(y_imp,
                       X_imp){

  # ----------------------------- preparing the X data ------------------
  # remove excessive variables
  X_imp <- remove_excessives(X_imp)

  # standardize X
  X_imp_stand <- stand(X_imp)

  #the missing indactor indicates, which values of y are missing.
  missind <- is.na(y_imp)

  n <- length(y_imp)

  #starting model
  ph <- sample_imp(y_imp)[, 1]
  tmp_0_all <- data.frame(target = ph)
  xnames_0 <- paste("X", 1:ncol(X_imp_stand), sep = "")
  tmp_0_all[xnames_0] <- X_imp_stand
  tmp_0_sub <- tmp_0_all[!missind, , drop = FALSE]

  reg_1_all <- stats::glm(target ~ 0 + ., data =  tmp_0_all, family = "poisson")
  reg_1_sub <- stats::glm(target ~ 0 + ., data =  tmp_0_sub, family = "poisson")

  X_model_matrix_1_all <- stats::model.matrix(reg_1_all)
  xnames_1 <- paste("X", 1:ncol(X_model_matrix_1_all), sep = "")

  #remove unneeded variables
  unneeded <- is.na(stats::coefficients(reg_1_sub))
  xnames_2 <- xnames_1[!unneeded]

  tmp_2_all <- data.frame(target = y_imp)
  X_model_matrix_2_all <- X_model_matrix_1_all[, !unneeded, drop = FALSE]

  tmp_2_all[, xnames_2] <- X_model_matrix_2_all

  tmp_2_sub <- tmp_2_all[!missind, , drop = FALSE]

  # -------------- calling the gibbs sampler to get imputation parameters----

  fixformula <- stats::formula(paste("target~", paste(xnames_2, collapse = "+"), "- 1",
                                     sep = ""))

  prior <- list(R = list(V = 1e-07, nu = -2))

  MCMCglmm_draws <- MCMCglmm::MCMCglmm(fixformula, data = tmp_2_all,
                                       verbose = FALSE, pr = TRUE, prior = prior,
                                       family = "poisson",
                                       saveX = TRUE,
                                       nitt = 3000,
                                       thin = 10,
                                       burnin = 1000)

  pointdraws <- MCMCglmm_draws$Sol
  xdraws <- pointdraws[, 1:ncol(X_model_matrix_2_all), drop = FALSE]
  variancedraws <- MCMCglmm_draws$VCV
  # the last column contains the variance (not standard deviation) of the residuals

  number_of_draws <- nrow(pointdraws)
  select.record <- sample(1:number_of_draws, 1, replace = TRUE)

  # -------------------- drawing samples with the parameters from the gibbs sampler --------
  ###start imputation

  fix.eff.imp <- matrix(xdraws[select.record, ], nrow = ncol(X_model_matrix_2_all))

  sigma.y.imp <- sqrt(variancedraws[select.record, ncol(variancedraws)])


  lambda <- exp(stats::rnorm(n, X_model_matrix_2_all %*% fix.eff.imp, sigma.y.imp))

  # note: the maximum lambda is something about 2.14e9
  if(max(lambda) > 2.14e9) {
    stop("Imputation of count variable failed due to a too high value of lambda.
This can occur if an observation in your data is an outlier regarding the covariates of the imputation model.
What again can be caused by highly variant imputation models for these covariates due high missing rates.")
  }
  y_proposal <- stats::rpois(n, lambda)
  y_tmp <- ifelse(missind, y_proposal, y_imp)
  y_ret <- data.frame(y_imp = y_tmp)


  # --------- returning the imputed data --------------
  return(y_ret)
}



# Generate documentation with devtools::document()
# Build package with devtools::build() and devtools::build(binary = TRUE) for zips
