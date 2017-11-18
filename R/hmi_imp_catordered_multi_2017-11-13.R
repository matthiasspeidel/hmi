#' The function for hierarchical imputation of categorical variables.
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
imp_orderedcat_multi <- function(y_imp,
                      X_imp,
                      Z_imp,
                      clID,
                      nitt = 25000,
                      burnin = 5000,
                      thin = 20){

  # ----- checking input
  if(!is.factor(y_imp)){
    warning("We suggest to make all your categorical variables to be a factor.")
    y_imp <- as.factor(y_imp)
  }

  # -----------------------------preparing the data ------------------
  X_imp <- cleanup(X_imp)
  X_imp_stand <- stand(X_imp)
  Z_imp <- cleanup(Z_imp)
  Z_imp_stand <- stand(Z_imp)

  ph <- sample_imp(y_imp)[, 1]

  missind <- is.na(y_imp)

  tmp_0_all <- data.frame(target = ph)
  xnames_0 <- paste("X", 1:ncol(X_imp_stand), sep = "")

  tmp_0_all[xnames_0] <- X_imp_stand
  tmp_0_sub <- tmp_0_all[!missind, , drop = FALSE]

  # -----------------------------gather some parameters
  n <- nrow(X_imp_stand)
  missind <- is.na(y_imp)
  n.par.rand <- ncol(Z_imp_stand)   # Get the number of random effects variables
  length.alpha <- length(table(clID)) * n.par.rand
  number_clusters <- length(table(clID))

  # ----------------- starting model to get categorical variables as dummies
  # and for eliminating linear dependent variables.
  reg_0_all <- stats::lm(as.numeric(target) ~ 0 + ., data =  tmp_0_all)
  reg_0_sub <- stats::lm(as.numeric(target) ~ 0 + ., data =  tmp_0_sub)

  unneeded_0 <- is.na(stats::coef(reg_0_sub))
  MM_0_all <- stats::model.matrix(reg_0_all)[, !unneeded_0, drop = FALSE]

  tmp_1_all <- data.frame(target = ph)
  xnames_1 <- paste("X", 1:ncol(MM_0_all), sep = "")
  tmp_1_all[xnames_1] <- MM_0_all
  tmp_1_sub <- tmp_1_all[!missind, , drop = FALSE]

  # -------- better model to remove insignificant variables

  types <- array(dim = ncol(tmp_0_all))
  for(i in 1:length(types)) types[i] <- get_type(tmp_0_all[, i])
  #polr needs an intercept variable. So we have to include one to the model formula,
  #but not if the data already have an intercept variable.

  tmp_1_all_nointercept <- tmp_1_all[, types != "intercept", drop = FALSE]
  tmp_1_sub_nointercept <- tmp_1_all[, types != "intercept", drop = FALSE]

  reg_1_all <- MASS::polr(target ~ 1 + ., data =  tmp_1_all_nointercept, method = "probit", Hess = TRUE)
  reg_1_sub <- MASS::polr(target ~ 1 + ., data =  tmp_1_sub_nointercept, method = "probit", Hess = TRUE)

  p <- ncol(stats::model.matrix(reg_1_all))
  insignificant <- c(FALSE, 2 * stats::pt(q = abs(summary(reg_1_sub)[[1]][1:(p-1), 3]),
                                   df = reg_1_sub$df.residual, lower.tail = FALSE) > 0.1)

  #---------------removing insignificant variables
  while(any(insignificant)){
    MM_1_all <- stats::model.matrix(reg_1_all)[, !insignificant, drop = FALSE]
    MM_1_sub <- MM_1_all[!missind, , drop = FALSE]

    tmp_2_all <- data.frame(target = ph)
    xnames_2 <- paste("X", 1:ncol(MM_1_all), sep = "")
    tmp_2_all[xnames_2] <- MM_1_all
    tmp_2_sub <- tmp_2_all[!missind, , drop = FALSE]

    types <- array(dim = ncol(tmp_2_all))
    for(i in 1:length(types)) types[i] <- get_type(tmp_2_all[, i])

    tmp_2_all_nointercept <- tmp_2_all[, types != "intercept", drop = FALSE]
    tmp_2_sub_nointercept <- tmp_2_all[, types != "intercept", drop = FALSE]

    reg_1_all <- MASS::polr(target ~ 1 + ., data =  tmp_2_all_nointercept, method = "probit", Hess = TRUE)
    reg_1_sub <- MASS::polr(target ~ 1 + ., data =  tmp_2_sub_nointercept, method = "probit", Hess = TRUE)

    p <- ncol(stats::model.matrix(reg_1_all))
    insignificant <- c(FALSE, 2 * stats::pt(q = abs(summary(reg_1_sub)[[1]][1:(p - 1), 3]),
                                     df = reg_1_sub$df.residual, lower.tail = FALSE) > 0.1)

  }


  MM_2_all <- stats::model.matrix(reg_1_all)
  xnames_2 <- paste("X", 1:ncol(MM_2_all), sep = "")
  znames <- paste("Z", 1:ncol(Z_imp_stand), sep = "")


  tmp_3_all <- data.frame(target = as.factor(y_imp))
  tmp_3_all[, xnames_2] <- MM_2_all
  tmp_3_all[, znames] <- Z_imp_stand
  tmp_3_all[, "clID"] <- clID
  tmp_3_sub <- tmp_3_all[!missind, , drop = FALSE]

  # -------------- calling the gibbs sampler to get imputation parameters----
  fixformula <- stats::formula(paste("target~ - 1 + ",
                                       paste(xnames_2, sep = "", collapse = "+")))

  randformula <- stats::formula(paste("~us( - 1 + ", paste(znames, sep = "", collapse = "+"),
                                 "):clID", sep = ""))

  number_fix_parameters <- ncol(MM_2_all)

  # Get the number of random effects variables
  number_random_effects <- length(znames)
  number_random_parameters <- number_random_effects
  #Fix residual variance R at 1
  # cf. http://stats.stackexchange.com/questions/32994/what-are-r-structure-g-structure-in-a-glmm

  J <- length(table(y_imp)) #number of categories
  #priors from https://stat.ethz.ch/pipermail/r-sig-mixed-models/2012q2/018194.html

  #The residual variance has to be fixed, because the latent variable is without scale
  prior <- list(R = list(V = 1, fix = TRUE),
       G = list(G1 = list(V = diag(number_random_effects), nu = 0.002)))

  MCMCglmm_draws <- MCMCglmm::MCMCglmm(fixed = fixformula,
                                       random = randformula,
                                       data = tmp_3_sub,
                                       family = "ordinal",
                                       verbose = FALSE, pr = TRUE,
                                       prior = prior,
                                       saveX = TRUE, saveZ = TRUE,
                                       nitt = nitt,
                                       thin = thin,
                                       burnin = burnin)

  pointdraws <- MCMCglmm_draws$Sol

  xdraws <- pointdraws[, 1:number_fix_parameters, drop = FALSE]
  zdraws <- pointdraws[, number_fix_parameters +
                         1:(number_random_parameters * number_clusters), drop = FALSE]

  variancedraws <- MCMCglmm_draws$VCV

  number_of_draws <- nrow(pointdraws)
  select.record <- sample(1:number_of_draws, size = 1)

  # -------------------- drawing samples with the parameters from the gibbs sampler --------

  #now decide in with category a person falls acording to threshold.

  ###start imputation
  rand.eff.imp <- matrix(zdraws[select.record,],
                           ncol = number_random_effects) # cluster specific effects

  fix.eff.imp <- matrix(xdraws[select.record, ], nrow = ncol(MM_2_all))

  sigma.y.imp <- sqrt(variancedraws[select.record, ncol(variancedraws)])

  latent <- stats::rnorm(n, MM_2_all %*% fix.eff.imp +
                      apply(Z_imp_stand * rand.eff.imp[clID,], 1, sum), sigma.y.imp)

    #cutpoint: the first cutpoint is set to 0 by MCMCglmm
    #(see https://stat.ethz.ch/pipermail/r-sig-mixed-models/2013q1/020030.html)

  cps <- c(min(latent), 0, MCMCglmm_draws$CP[select.record,], max(latent))

  ret0 <- y_imp
  ret0[missind] <- levels(y_imp)[as.numeric(cut(latent, breaks = cps, include.lowest = TRUE))][missind]
  y_ret <- data.frame(y_ret = ret0)

  # --------- returning the imputed data --------------
  ret <- list(y_ret = y_ret, Sol = xdraws, VCV = variancedraws)
  return(ret)
}
