#' The function for hierarchical imputation of continuous variables.
#'
#' The function is called by the wrapper.
#' @param y_imp A vector with the variable to impute.
#' @param X_imp A data.frame with the fixed effects variables.
#' @param Z_imp A data.frame with the random effects variables.
#' @param clID A vector with the cluster ID.
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
imp_cont_multi <- function(y_imp,
                      X_imp,
                      Z_imp,
                      clID,
                      nitt = 22000,
                      burnin = 2000,
                      thin = 20,
                      pvalue = 0.2,
                      rounding_degrees = c(1, 10, 100, 1000)){

  # The missing indactor indicates, which values of y are missing.
  missind <- is.na(y_imp)
  n <- length(y_imp)

  # -----------------------------preparing the data ------------------
  # -- standardise the covariates in X (which are numeric and no intercept)
  # ----------------------------- preparing the X and Z data ------------------

  # remove excessive variables
  X_imp <- cleanup(X_imp)

  # standardise the covariates in X (which are numeric and no intercept)
  X_imp_stand <- stand(X_imp, rounding_degrees = rounding_degrees)

  # -- standardise the covariates in Z (which are numeric and no intercept)
  Z_imp <- cleanup(Z_imp)

  Z <- stand(Z_imp, rounding_degrees = rounding_degrees)

  #define a place holder (ph)
  ph <- sample_imp(y_imp)[, 1]

  y_mean <- mean(ph, na.rm = TRUE)
  y_sd <- stats::sd(ph, na.rm = TRUE)

  ph <- (ph - y_mean)/y_sd + 1

  YZ <- data.frame(target = ph, Z)
  #remove intercept variable
  YZ <- YZ[, apply(YZ, 2, get_type) != "intercept", drop = FALSE]

  Z2 <- stats::model.matrix(stats::lm("target ~ 1 + .", data = YZ))

  # ----------- set up a maximal model matrix with all possible relevant (dummy) variables -----
  # In the imputation model only actually relevant (dummy) variables shall be present.
  # THis is done by setting up a mirror of the initial model matrix.
  # Then step by step this model matrix is reduced to all actually relevant (dummy) variables.
  # This reduction is based on models using the observed data.
  # The last step prior to the imputation-parameters estimation is to restrict the initial mode matrix
  # to those variables, left in the reduced mirror model matrix.


  tmp_0_all <- data.frame(target = ph)
  xnames_1 <- paste("X", 1:ncol(X_imp_stand), sep = "")
  znames_1 <- paste("Z", 1:ncol(Z2), sep = "")

  tmp_0_all[, xnames_1] <- X_imp_stand
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
  suppressMessages(reg_1_all <- lme4::lmer(stats::formula(tmp_formula), data = tmp_0_all))


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
  reg_1_sub <- lme4::lmer(stats::formula(tmp_formula) , data = tmp_1_sub)

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

    if(length(xnames_1) == 0){
      fixformula_lme <- stats::formula("target~ 1")
    }else{
      fixformula_lme <- stats::formula(paste("target~ 0+ ", paste(xnames_1, collapse = "+"), sep = ""))
    }

    if(length(znames_1) == 0){
      randformula_lme <- stats::as.formula("~1|clID")
    }else{
      randformula_lme <- stats::as.formula(paste("~0+", paste(znames_1, collapse = "+"), "|clID",
                                                 sep = ""))
    }

    reg_1_sub <- nlme::lme(fixed = fixformula_lme,
                           random = randformula_lme, data = tmp_1_sub)

    pvalues <- stats::anova(reg_1_sub)$"p-value"
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

  # -------------- calling the gibbs sampler to get imputation parameters----

  prior <- list(R = list(V = 1, nu = 0.002), # alternative: R = list(V = 1e-07, nu = -2)
                G = list(G1 = list(V = diag(ncol(Z2)), nu = 0.002)))

  #run MCMCglmm based on the data with observations not missing and variables not unimportant
  MCMCglmm_draws <- MCMCglmm::MCMCglmm(fixed = fixformula, random = randformula,
                                       data = YXZ_2_sub,
                                       verbose = FALSE, pr = TRUE, prior = prior,
                                       saveX = TRUE, saveZ = TRUE,
                                       nitt = nitt,
                                       thin = thin,
                                       burnin = burnin)

  tmp_2_all <- tmp_0_all[, colnames(YXZ_2_sub), drop = FALSE]
  # Get the number of random effects variables
  n.par.rand <- ncol(Z2)
  ncluster <- length(table(YXZ_2_sub$clID))
  length.alpha <- ncluster * n.par.rand

  pointdraws <- MCMCglmm_draws$Sol
  xdraws <- pointdraws[, 1:ncol(X_model_matrix_1_sub), drop = FALSE]
  #If a cluster cannot has random effects estimates, because there too few observations,
  #we make them 0.
  empty_cluster <- which(table(YXZ_2_sub$clID) == 0)
  zdraws_pre <- pointdraws[, (ncol(X_model_matrix_1_sub) + 1):ncol(pointdraws), drop = FALSE]
  #go through all random effects
  #(e.g. first the random intercepts, then the random slope of X1, then the random slope of X5)
  for(l1 in 1:n.par.rand){
    #go through all clusters with 0 observations
    for(l2 in empty_cluster){
      zdraws_pre <- cbind(zdraws_pre[, 0:((l1-1)* ncluster + (l2-1))], 0,
                          zdraws_pre[,  ((l1-1)* ncluster + l2):ncol(zdraws_pre)])
    }
  }

  zdraws <- zdraws_pre
  variancedraws <- MCMCglmm_draws$VCV
  # the last column contains the variance (not standard deviation) of the residuals

  number_of_draws <- nrow(pointdraws)
  select.record <- sample(1:number_of_draws, size = 1)

  # -------------------- drawing samples with the parameters from the gibbs sampler --------
  ###start imputation
  rand.eff.imp <- matrix(zdraws[select.record, ], ncol = n.par.rand)

  fix.eff.imp <- matrix(xdraws[select.record, ], nrow = ncol(X_model_matrix_1_sub))

  sigma.y.imp <- sqrt(variancedraws[select.record, ncol(variancedraws)])

  y_temp <- stats::rnorm(n, as.matrix(tmp_2_all[, xnames_1, drop = FALSE])  %*% fix.eff.imp +
                      apply(Z2 * rand.eff.imp[clID, , drop = FALSE], 1, sum), sd = sigma.y.imp)

  y_ret <- data.frame(y_ret = ifelse(missind, (y_temp - 1) * y_sd + y_mean, y_imp))

  # --------- returning the imputed data --------------
  ret <- list(y_ret = y_ret, Sol = xdraws, VCV = variancedraws)
  return(ret)
}
