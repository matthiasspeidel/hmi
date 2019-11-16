#' The function for hierarchical imputation of binary variables.
#'
#' The function is called by the wrapper.
#' @param y_imp A Vector with the variable to impute.
#' @param X_imp A data.frame with the fixed effects variables.
#' @param Z_imp A data.frame with the random effects variables.
#' @param clID A vector with the cluster ID.
#' @param nitt An integer defining number of MCMC iterations (see MCMCglmm).
#' @param burnin burnin A numeric value between 0 and 1 for the desired percentage of
#' Gibbs samples that shall be regarded as burnin.
#' @param thin An integer to set the thinning interval range. If thin = 1,
#' every iteration of the Gibbs-sampling chain will be kept. For highly autocorrelated
#' chains, that are only examined by few iterations (say less than 1000),
#' @param pvalue A numeric between 0 and 1 denoting the threshold of p-values a variable in the imputation
#' model should not exceed. If they do, they are excluded from the imputation model.
#' @param k An integer defining the allowed maximum of levels in a factor covariate.
#' @return A list with 1. 'y_ret' the n x 1 data.frame with the original and imputed values.
#' 2. 'Sol' the Gibbs-samples for the fixed effects parameters.
#' 3. 'VCV' the Gibbs-samples for variance parameters.
imp_binary_multi <- function(y_imp,
                      X_imp,
                      Z_imp,
                      clID,
                      nitt = 22000,
                      burnin = 2000,
                      thin = 20,
                      pvalue = 0.2,
                      k = Inf){

  # ----------------------------- preparing the y data ------------------
  # stransform y_imp into a real binary with only zeros and ones (and NAs).
  first_possibility <- sort(y_imp)[1]
  second_possibility <- rev(sort(y_imp))[1]
  y_binary <- data.frame(y = factor(y_imp, labels = c(0, 1)))


  # If one category has less then two observations, no binary model can be estimated.
  # So the imputation routines has to stop.
  if(min(table(y_binary)) < 2){
    stop("A binary (or maybe a semicontinuous) variable has less than two observations in one category.
         Consider removing this variable
         (or in the case of a semicontinuous variable, to specify it as continuous in the list_of_types (see ?hmi)).")
  }

  # The missing indactor indicates, which values of y are missing.
  missind <- is.na(y_binary)
  n <- length(y_imp)

  # ----------------------------- preparing the X and Z data ------------------

  # remove excessive variables
  X_imp <- cleanup(X_imp, k = k)

  # standardise the covariates in X (which are numeric and no intercept)
  X <- stand(X_imp)

  # -- standardise the covariates in Z (which are numeric and no intercept)
  Z_imp <- cleanup(Z_imp, k = k)
  Z <- stand(Z_imp)


  # Get the number of random effects variables
  n.par.rand <- ncol(Z)
  length.alpha <- length(table(clID)) * n.par.rand

  # We need two design matrices. One for the model to get the imputation parameters
  # which is based on the observed values only (obs).
  # And secondly a design matrix for E(y_mis|X_mis) based on all observations (all).
  # The later is more or less only the technical framework to get the imputed values.

  #define a place holder (ph)
  ph <- sample_imp(y_binary[, 1])[, 1]

  YZ <- data.frame(target = ph, Z)
  #remove intercept variable
  YZ <- YZ[, apply(YZ, 2, get_type) != "intercept", drop = FALSE]

  Z2 <- stats::model.matrix(stats::glm("target ~ 1 + .", data = YZ,
                                       family = stats::binomial(link = "logit")))
  # ----------- set up a maximal model matrix with all possible relevant (dummy) variables -----
  # In the imputation model only actually relevant (dummy) variables shall be present.
  # THis is done by setting up a mirror of the initial model matrix.
  # Then step by step this model matrix is reduced to all actually relevant (dummy) variables.
  # This reduction is based on models using the observed data.
  # The last step prior to the imputation-parameters estimation is to restrict the initial mode matrix
  # to those variables, left in the reduced mirror model matrix.


  tmp_0_all <- data.frame(target = ph)
  xnames_1 <- paste("X", 1:ncol(X), sep = "")
  znames_1 <- paste("Z", 1:ncol(Z2), sep = "")

  tmp_0_all[, xnames_1] <- X
  tmp_0_all[, znames_1] <- Z2
  tmp_0_all[, "clID"] <- clID

  tmp_formula <- paste("target ~ 0 +",
                       paste(xnames_1, collapse = "+"),
                       "+(0+",
                       paste(znames_1, collapse = "+"),
                       "|clID)")

  # If both, an intercept variable and a categorical variable are present in the data,
  # One variable in the model is redundant. This is handled later in the code, so here
  # the default message from lmer is bothering and therefore suppressed.
  oldw <- getOption("warn")
  options(warn = -1)
  suppressMessages(reg_1_all <- lme4::glmer(stats::formula(tmp_formula),
                                            family = stats::binomial("logit"), data = tmp_0_all))
  options(warn = oldw)

  X_model_matrix_1_all <- stats::model.matrix(reg_1_all)
  xnames_1 <- paste("X", 1:ncol(X_model_matrix_1_all), sep = "")
  colnames(X_model_matrix_1_all) <- xnames_1

  tmp_0_all <- data.frame(target = ph)
  tmp_0_all[, xnames_1] <- X_model_matrix_1_all
  tmp_0_all[, znames_1] <- Z2
  tmp_0_all[, "clID"] <- clID

  #From this initial model matrix X_model_matrix_1_all
  #now step by step irrelavant variables are removed.
  X_model_matrix_1_sub <- X_model_matrix_1_all[!missind, , drop = FALSE]

  # The first step of the reduction is to remove variables having a non-measurable effect
  # (e.g. due to colinearity) on y.
  # tmp_1 shall include the covariates (like X_model_matrix) and additionally the target variable
  ph_sub <- ph[!missind]
  tmp_1_sub <- data.frame(target = ph_sub)
  tmp_1_sub[, xnames_1] <- X_model_matrix_1_sub
  tmp_1_sub[, znames_1] <- Z2[!missind, , drop = FALSE]
  tmp_1_sub[, "clID"] <- clID[!missind]

  tmp_formula <- paste("target ~ 0 +",
                       paste(xnames_1, collapse = "+"),
                       "+(0+",
                       paste(znames_1, collapse = "+"),
                       "|clID)")

  oldw <- getOption("warn")
  options(warn = -1)
  suppressMessages(reg_1_sub <- lme4::glmer(stats::formula(tmp_formula),
                                            family = stats::binomial("logit"), data = tmp_1_sub))
  options(warn = oldw)

  #remove unneeded variables
  X_model_matrix_1_sub <- X_model_matrix_1_sub[, !is.na(lme4::fixef(reg_1_sub)),
                                               drop = FALSE]

  # Remove insignificant variables from the imputation model
  check <- TRUE
  while(check){
    tmp_1_sub <- data.frame(target = ph_sub)
    xnames_1 <- colnames(X_model_matrix_1_sub)
    tmp_1_sub[, xnames_1] <- X_model_matrix_1_sub
    tmp_1_sub[, znames_1] <- Z2[!missind, , drop = FALSE]
    tmp_1_sub[, "clID"] <- clID[!missind]

    tmp_formula <- paste("target ~ 0 +",
                         paste(xnames_1, collapse = "+"),
                         "+(0+",
                         paste(znames_1, collapse = "+"),
                         "|clID)")

    oldw <- getOption("warn")
    options(warn = -1)
    suppressMessages(reg_1_sub <- lme4::glmer(stats::formula(tmp_formula),
                                              family = stats::binomial("logit"), data = tmp_1_sub))
    options(warn = oldw)

    pvalues <- summary(reg_1_sub)$coefficients[, "Pr(>|z|)"]
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


  YXZ_2_sub <- data.frame(target = ph_sub)
  xnames_1 <- colnames(X_model_matrix_1_sub)
  YXZ_2_sub[, xnames_1] <- X_model_matrix_1_sub
  YXZ_2_sub[, znames_1] <- Z2[!missind, , drop = FALSE]
  YXZ_2_sub[, "clID"] <- clID[!missind]

  if(length(xnames_1) == 0){
    fixformula <- stats::formula("target~ 1")
  }else{
    fixformula <- stats::formula(paste("target~ 0 + ", paste(xnames_1, collapse = "+"), sep = ""))
  }

  if(length(znames_1) == 0){
    randformula <- stats::as.formula("~us(1):ID")
  }else{
    randformula <- stats::as.formula(paste("~us(0 +", paste(znames_1, collapse = "+"), "):clID",
                                           sep = ""))
  }

  # -------------- calling the gibbs sampler to get imputation parameters ----

  #Fix residual variance R at 1
  # cf. http://stats.stackexchange.com/questions/32994/what-are-r-structure-g-structure-in-a-glmm
  prior <- list(R = list(V = 1, fix = TRUE),
                G = list(G1 = list(V = diag(n.par.rand), nu = 0.002)))

  MCMCglmm_draws <- MCMCglmm::MCMCglmm(fixed = fixformula,
                                       random = randformula,
                                       data = YXZ_2_sub,
                                       family = "categorical",
                                       verbose = FALSE, pr = TRUE, prior = prior,
                                       saveX = TRUE, saveZ = TRUE,
                                       nitt = nitt,
                                       thin = thin,
                                       burnin = burnin)

  tmp_2_all <- tmp_0_all[, colnames(YXZ_2_sub), drop = FALSE]
  # correction. see:
  # http://stats.stackexchange.com/questions/32994/what-are-r-structure-g-structure-in-a-glmm
  k <- ((16*sqrt(3))/(15*pi))^2

  pointdraws <- MCMCglmm_draws$Sol /sqrt(1 + k)
  xdraws <- pointdraws[, 1:ncol(X_model_matrix_1_sub), drop = FALSE]
  zdraws <- pointdraws[, ncol(X_model_matrix_1_sub) + 1:length.alpha, drop = FALSE]
  variancedraws <- MCMCglmm_draws$VCV

  number_of_draws <- nrow(pointdraws)
  select_record <- sample(1:number_of_draws, size = 1)

  # -------------------- drawing samples with the parameters from the gibbs sampler --------

  linkfunction <- boot::inv.logit

  ###start imputation

  rand_eff_imp <- matrix(zdraws[select_record,],
                           ncol = n.par.rand)

  fix_eff_imp <- matrix(xdraws[select_record, ], nrow = ncol(X_model_matrix_1_sub))

  sigma_y_imp <- sqrt(variancedraws[select_record, ncol(variancedraws)])

  linearpredictor <- stats::rnorm(n, as.matrix(tmp_2_all[, xnames_1, drop = FALSE]) %*% fix_eff_imp +
                      apply(Z2 * rand_eff_imp[clID,], 1, sum), 0*sigma_y_imp)

  #chances to draw a get a 1 (which means to get the second_possibility)
  one_prob <- linkfunction(linearpredictor)

  # determine whether an observation gets the 1 according the the chance
  # calculated in one_prob.
  indicator <- as.numeric(stats::runif(n) < one_prob)

  #Initialising the returning matrix

  y_ret <- as.data.frame(y_imp)
  for(i in which(is.na(y_imp))){
    if(indicator[i] == 0){
      y_ret[i, 1] <- first_possibility
    }else{
      y_ret[i, 1] <- second_possibility
    }
  }
  colnames(y_ret) <- "y_ret"

  # --------- returning the imputed data --------------
  ret <- list(y_ret = y_ret, Sol = xdraws, VCV = variancedraws)
  return(ret)
}
