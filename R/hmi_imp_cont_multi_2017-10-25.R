#' The function for hierarchical imputation of continuous variables.
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
#' @return A list with 1. 'y_ret' the n x 1 data.frame with the original and imputed values.
#' 2. 'Sol' the Gibbs-samples for the fixed effects parameters.
#' 3. 'VCV' the Gibbs-samples for variance parameters.
imp_cont_multi <- function(y_imp,
                      X_imp,
                      Z_imp,
                      clID,
                      nitt = 22000,
                      burnin = 2000,
                      thin = 20){

  # -----------------------------preparing the data ------------------
  # -- standardise the covariates in X (which are numeric and no intercept)
  # ----------------------------- preparing the X and Z data ------------------

  # remove excessive variables
  X_imp <- cleanup(X_imp)

  # standardise the covariates in X (which are numeric and no intercept)
  X <- stand(X_imp)

  # -- standardise the covariates in Z (which are numeric and no intercept)
  Z_imp <- cleanup(Z_imp)

  Z <- stand(Z_imp)


  #define a place holder (ph)
  ph <- sample_imp(y_imp)[, 1]

  y_mean <- mean(ph, na.rm = TRUE)
  y_sd <- stats::sd(ph, na.rm = TRUE)

  ph_stand <- (ph - y_mean)/y_sd + 1
  #Z may contain factor variables. Later we need Z to contain numeric variables,
  # so we have to change those variables.
  YZ <- data.frame(target = ph_stand, Z)
  #remove intercept variable
  YZ <- YZ[, get_type(YZ) != "intercept", drop = FALSE]

  Z2 <- stats::model.matrix(stats::lm("target ~ 1 + .", data = YZ))

  missind <- is.na(y_imp)
  n <- nrow(X)

  # ----------------- starting model to get categorical variables as dummies
  # and for eliminating linear dependent variables.
  #get data set with all variables and all observations
  xnames_0 <- paste("X", 1:ncol(X), sep = "")
  YX_0_all <- data.frame(target = ph_stand)
  YX_0_all[xnames_0] <- X
  YX_0_sub <- YX_0_all[!missind, , drop = FALSE]

  reg_YX_0_all <- stats::lm(target ~ 0 + ., data =  YX_0_all)
  reg_YX_0_sub <- stats::lm(target ~ 0 + ., data =  YX_0_sub)

  unneeded_0 <- is.na(stats::coef(reg_YX_0_sub))
  X_0_all <- stats::model.matrix(reg_YX_0_all)[, !unneeded_0, drop = FALSE]

  YXZ_1_all <- data.frame(target = ph_stand)
  xnames_1 <- paste("X", 1:ncol(X_0_all), sep = "")
  znames_1 <- paste("Z", 1:ncol(Z2), sep = "")
  YXZ_1_all[, xnames_1] <- X_0_all
  YXZ_1_all[, znames_1] <- Z2
  YXZ_1_all[, "clID"] <- clID
  YXZ_1_sub <- YXZ_1_all[!missind, , drop = FALSE]

  # -------- better model to remove insignificant variables
  formula_0 <- stats::as.formula(paste("target~ 0 +",
                                paste(xnames_1, collapse = "+"),
                                "+(0+",
                                paste(znames_1, collapse = "+"),
                                "|clID)"))

  reg_YXZ_1_all <- lme4::lmer(formula_0, data = YXZ_1_all)
  reg_YXZ_1_sub <- lme4::lmer(formula_0, data = YXZ_1_sub)
  #note that lmer does not provide p-values or degrees of freedom.
  #We follow the approach to accept variables that are insignificant in the imputation model,
  #so higher degrees of freedom yield more results "parameter significantly different from 0".
  #As a rule of thumb, 30 degrees of freedom approximate the normal distribution quite well
  #and for n -> Inf, the t-distribution approaches the normal distribution.
  #so we will use the normal distribution to find variables that seem to be not significant.
  insignificant <- 2 * stats::pnorm(q = abs(summary(reg_YXZ_1_sub)[[10]][, 3]), lower.tail = FALSE) > 0.1

  while(any(insignificant)){
    X_0_all <- stats::model.matrix(reg_YXZ_1_all)[, !insignificant, drop = FALSE]

    YXZ_1_all <- data.frame(target = ph_stand)
    xnames_1 <- paste("X", 1:ncol(X_0_all), sep = "")
    YXZ_1_all[xnames_1] <- X_0_all
    YXZ_1_all[, znames_1] <- Z2
    YXZ_1_all[, "clID"] <- clID
    YXZ_1_sub <- YXZ_1_all[!missind, , drop = FALSE]

    formula_0 <- stats::as.formula(paste("target ~ 0 +",
                                  paste(xnames_1, collapse = "+"),
                                  "+(0+",
                                  paste(znames_1, collapse = "+"),
                                  "|clID)"))

    reg_YXZ_1_all <- lme4::lmer(formula_0, data = YXZ_1_all)
    reg_YXZ_1_sub <- lme4::lmer(formula_0, data = YXZ_1_sub)

    insignificant <- 2 * stats::pnorm(q = abs(summary(reg_YXZ_1_sub)[[10]][, 3]), lower.tail = FALSE) > 0.1
  }


  #drop intercept variables from the data.set
  intercept_variables <- sapply(YXZ_1_all, get_type) == "intercept"
  YXZ_2_all <- YXZ_1_all[, !intercept_variables, drop = FALSE]
  YXZ_2_sub <- YXZ_2_all[!missind, , drop = FALSE]

  #ensure, that intercept variables do not appear in the fixformula or randformula

  xnames_2 <- xnames_1[!xnames_1 %in% names(intercept_variables)[intercept_variables]]
  znames_2 <- znames_1[!znames_1 %in% names(intercept_variables)[intercept_variables]]

  if(length(xnames_2) == 0){
    fixformula <- stats::formula("target~ 1")
  }else{
    fixformula <- stats::formula(paste("target~ 1+ ", paste(xnames_2, collapse = "+"), sep = ""))
  }

  if(length(znames_2) == 0){
    randformula <- stats::as.formula("~us(1):ID")
  }else{
    randformula <- stats::as.formula(paste("~us(1+", paste(znames_2, collapse = "+"), "):clID",
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

  # Get the number of random effects variables
  n.par.rand <- ncol(Z2)
  ncluster <- length(table(YXZ_2_sub$clID))
  length.alpha <- ncluster * n.par.rand

  pointdraws <- MCMCglmm_draws$Sol
  xdraws <- pointdraws[, 1:ncol(X_0_all), drop = FALSE]
  #If a cluster cannot has random effects estimates, because there too few observations,
  #we make them 0.
  empty_cluster <- which(table(YXZ_2_all$clID) == 0)
  zdraws_pre <- pointdraws[, (ncol(X_0_all) + 1):ncol(pointdraws), drop = FALSE]
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

  fix.eff.imp <- matrix(xdraws[select.record, ], nrow = ncol(X_0_all))

  sigma.y.imp <- sqrt(variancedraws[select.record, ncol(variancedraws)])

  y_temp <- stats::rnorm(n, X_0_all %*% fix.eff.imp +
                      apply(Z2 * rand.eff.imp[clID, ], 1, sum), sd = sigma.y.imp)

  y_ret <- data.frame(y_ret = ifelse(missind, (y_temp - 1) * y_sd + y_mean, y_imp))

  # --------- returning the imputed data --------------
  ret <- list(y_ret = y_ret, Sol = xdraws, VCV = variancedraws)
  return(ret)
}
