#' The function for imputation of binary variables.
#'
#' The function is called by the wrapper.
#' @param y_imp A vector with the variable to impute.
#' @param X_imp A data.frame with the fixed effects variables.
#' @param nitt An integer defining number of MCMC iterations (see MCMCglmm).
#' @param burnin burnin A numeric value between 0 and 1 for the desired percentage of
#' Gibbs samples that shall be regarded as burnin.
#' @param thin An integer to set the thinning interval range. If thin = 1,
#' every iteration of the Gibbs-sampling chain will be kept. For highly autocorrelated
#' chains, that are only examined by few iterations (say less than 1000).
#' @param pvalue A numeric between 0 and 1 denoting the threshold of p-values a variable in the imputation
#' model should not exceed. If they do, they are excluded from the imputation model.
#' @param rounding_degrees A numeric vector with the presumed rounding degrees.
#' @return A list with 1. 'y_ret' the n x 1 data.frame with the original and imputed values.
#' 2. 'Sol' the Gibbs-samples for the fixed effects parameters.
#' 3. 'VCV' the Gibbs-samples for variance parameters.
imp_count_single <- function(y_imp,
                       X_imp,
                       nitt = 22000,
                       burnin = 2000,
                       thin = 20,
                       pvalue = 0.2,
                       rounding_degrees = c(1, 10, 100, 1000)){

  # ----------------------------- preparing the X data ------------------
  # remove excessive variables
  X_imp <- cleanup(X_imp)

  # standardize X
  X <- stand(X_imp, rounding_degrees = rounding_degrees)

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

  tmp_0_all <- data.frame(target = ph, X)
  xnames_1 <- colnames(X)

  tmp_formula <- paste("target~ 0 + ", paste(xnames_1, collapse = "+"), sep = "")
  reg_1_all <- stats::glm(stats::formula(tmp_formula), data = tmp_0_all,
                          family = "poisson")

  X_model_matrix_1_all <- stats::model.matrix(reg_1_all)
  xnames_1 <- paste("X", 1:ncol(X_model_matrix_1_all), sep = "")
  colnames(X_model_matrix_1_all) <- xnames_1

  tmp_0_all <- data.frame(target = ph)
  tmp_0_all[, xnames_1] <- X_model_matrix_1_all

  #From this initial model matrix X_model_matrix_1_all
  #now step by step irrelavant variables are removed.
  X_model_matrix_1_sub <- X_model_matrix_1_all[!missind, , drop = FALSE]


  #first step of the reduction is to remove variables having a NA-effect (e.g. due to colinearity) on y
  #tmp_1 shall include the covariates (like X_model_matrix) and additionally the target variable
  ph_sub <- y_imp[!missind]
  tmp_1_sub <- data.frame(target = ph_sub)
  xnames_1 <- colnames(X_model_matrix_1_sub)
  tmp_1_sub[, xnames_1] <- X_model_matrix_1_sub

  tmp_formula <- paste("target~ 0 + ", paste(xnames_1, collapse = "+"), sep = "")

  reg_1_sub <- stats::glm(stats::formula(tmp_formula), data = tmp_1_sub,
                          family = "poisson")

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
                            family = "poisson")

    pvalues <- summary(reg_1_sub)$coefficients[, 4]
    insignificant_variables <- which(pvalues > pvalue)
    most_insignificant <- insignificant_variables[which.max(pvalues[insignificant_variables])]

    if(length(most_insignificant) == 0){
      check <- FALSE
    }else{
      X_model_matrix_1_sub <- stats::model.matrix(reg_1_sub)[, -most_insignificant, drop = FALSE]
    }
  }

  tmp_2_all <- tmp_0_all[, colnames(tmp_1_sub), drop = FALSE]
  tmp_2_all$target[missind] <- NA

  # -------------- calling the gibbs sampler to get imputation parameters----

  fixformula <- stats::formula(paste("target~", paste(xnames_1, collapse = "+"), "- 1",
                                     sep = ""))

  prior <- list(R = list(V = 1e-07, nu = -2))

  MCMCglmm_draws <- MCMCglmm::MCMCglmm(fixformula, data = tmp_1_sub,
                                       verbose = FALSE, pr = TRUE, prior = prior,
                                       family = "poisson",
                                       saveX = TRUE,
                                       nitt = nitt,
                                       thin = thin,
                                       burnin = burnin)

  pointdraws <- MCMCglmm_draws$Sol
  xdraws <- pointdraws[, 1:ncol(X_model_matrix_1_sub), drop = FALSE]
  variancedraws <- MCMCglmm_draws$VCV
  # the last column contains the variance (not standard deviation) of the residuals

  number_of_draws <- nrow(pointdraws)
  select.record <- sample(1:number_of_draws, size = 1)

  # -------------------- drawing samples with the parameters from the gibbs sampler --------
  ###start imputation

  fix.eff.imp <- matrix(xdraws[select.record, ], nrow = ncol(X_model_matrix_1_sub))

  sigma.y.imp <- sqrt(variancedraws[select.record, ncol(variancedraws)])


  lambda <- exp(stats::rnorm(n, as.matrix(tmp_2_all[, xnames_1, drop = FALSE]) %*% fix.eff.imp, sigma.y.imp))

  # note: the maximum lambda is something about 2.14e9
  if(max(lambda) > 2.14e9) {
    stop("Imputation of count variable failed due to a too high value of lambda.
This can occur if an observation in your data is an outlier regarding the covariates of the imputation model.
What again can be caused by highly variant imputation models for these covariates due high missing rates.")
  }
  y_proposal <- stats::rpois(n, lambda)
  y_tmp <- ifelse(missind, y_proposal, y_imp)
  y_ret <- data.frame(y_ret = y_tmp)


  # --------- returning the imputed data --------------
  ret <- list(y_ret = y_ret, Sol = xdraws, VCV = variancedraws)
  return(ret)
}



# Generate documentation with devtools::document()
# Build package with devtools::build() and devtools::build(binary = TRUE) for zips
