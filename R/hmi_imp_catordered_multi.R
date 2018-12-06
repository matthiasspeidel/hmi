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
#' chains, that are only examined by few iterations (say less than 1000).
#' @param pvalue A numeric between 0 and 1 denoting the threshold of p-values a variable in the imputation
#' model should not exceed. If they do, they are excluded from the imputation model.
#' @param k An integer defining the allowed maximum of levels in a factor covariate.
#' @return A list with 1. 'y_ret' the n x 1 data.frame with the original and imputed values.
#' 2. 'Sol' the Gibbs-samples for the fixed effects parameters.
#' 3. 'VCV' the Gibbs-samples for variance parameters.
imp_orderedcat_multi <- function(y_imp,
                      X_imp,
                      Z_imp,
                      clID,
                      nitt = 25000,
                      burnin = 5000,
                      thin = 20,
                      pvalue = 0.2,
                      k = Inf){

  # ----- checking input
  if(!is.factor(y_imp)){
    warning("We suggest to make all your categorical variables to be a factor.")
    y_imp <- as.factor(y_imp)
  }


  # remove excessive variables
  X_imp <- cleanup(X_imp, k = k)
  # standardise the covariates in X (which are numeric and no intercept)
  X <- stand(X_imp)

  # -- standardise the covariates in Z (which are numeric and no intercept)
  Z_imp <- cleanup(Z_imp, k = k)
  Z <- stand(Z_imp)

  #the missing indactor indicates, which values of y are missing.
  missind <- is.na(y_imp)
  n <- nrow(X)
  number_clusters <- length(table(clID))


  # ----------- set up a maximal model matrix with all possible relevant (dummy) variables -----
  # In the imputation model only actually relevant (dummy) variables shall be present.
  # THis is done by setting up a mirror of the initial model matrix.
  # Then step by step this model matrix is reduced to all actually relevant (dummy) variables.
  # This reduction is based on models using the observed data.
  # The last step prior to the imputation-parameters estimation is to restrict the initial mode matrix
  # to those variables, left in the reduced mirror model matrix.
  ph <- sample_imp(y_imp)[, 1]

  YZ <- data.frame(target = rnorm(n), Z)
  #remove intercept variable
  YZ <- YZ[, apply(YZ, 2, get_type) != "intercept", drop = FALSE]

  Z2 <- stats::model.matrix(stats::lm("target ~ 1 + .", data = YZ))

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
  suppressMessages(reg_1_all <- ordinal::clmm(stats::formula(tmp_formula), data = tmp_0_all,
                                              model = TRUE, Hess = TRUE))
  options(warn = oldw)

  X_model_matrix_1_all <- stats::model.matrix(stats::formula(paste("target ~ 0 +",
                                paste(xnames_1, collapse = "+"))), data = tmp_0_all)
    #stats::model.matrix(reg_1_all)
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
  suppressMessages(reg_1_sub <- ordinal::clmm(stats::formula(tmp_formula), data = tmp_1_sub))
  options(warn = oldw)

  #remove unneeded variables
  tmp <- stats::coefficients(reg_1_sub)
  #clmm does not estimate one effect for the intercept variable,
  #but K - 1 for the thresholds.
  # So only variables in tmp, that are also found in X_model_matrix_1_sub are evaluated.
  if(length(names(tmp)[is.na(tmp)]) > 0){
    X_model_matrix_1_sub <- X_model_matrix_1_sub[, !names(X_model_matrix_1_sub) %in% names(tmp)[is.na(tmp)],
                                                 drop = FALSE]
  }


  ############################################################
  # Remove insignificant variables from the imputation model #
  ############################################################
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
    suppressMessages(reg_1_sub <- ordinal::clmm(stats::formula(tmp_formula), data = tmp_1_sub))
    options(warn = oldw)
    tmp <- summary(reg_1_sub)$coefficients[, "Pr(>|z|)"]
    pvalues <- tmp[names(tmp) %in% xnames_1]
    insignificant_variables <- which(pvalues > pvalue)
    most_insignificant <- insignificant_variables[which.max(pvalues[insignificant_variables])]

    if(length(most_insignificant) == 0){
      check <- FALSE
    }else{
      tmp <- stats::model.matrix(reg_1_sub) #if an additional intercept variable is included by the model
      #we cannot run stats::model.matrix(reg_1_sub)[, -most_insignificant]
      #Because most_insignificant refers to a situation without an intercept variable.

      tmp_MM <- tmp[, !colnames(tmp) %in% names(most_insignificant), drop = FALSE]
      if(ncol(tmp_MM) == 0){
        check <- FALSE
      }else{
        X_model_matrix_1_sub <- tmp_MM
      }
    }
  }

  #Remove highly correlated variables
  #Less save alternative: get_type(X_model_matrix_1_sub)
  tmptypes <- array(dim = ncol(X_model_matrix_1_sub))
  for(j in 1:ncol(X_model_matrix_1_sub)){

    tmptypes[j] <- get_type(X_model_matrix_1_sub[, j])
  }
  somethingcont <- tmptypes %in% c("cont", "binary", "roundedcont", "semicont")
  correlated <- NULL
  if(sum(somethingcont, na.rm = TRUE) >= 2){
    tmp0 <- stats::cor(X_model_matrix_1_sub[, somethingcont])
    tmp0[lower.tri(tmp0, diag = TRUE)] <- NA
    tmp1 <- stats::na.omit(as.data.frame(as.table(tmp0)))
    correlated <- unique(as.character(subset(tmp1,
                                             abs(tmp1[, 3]) > 0.99 & abs(tmp1[, 3]) < 1)[, 1]))
  }

  if(!is.null(correlated)){
    X_model_matrix_1_sub <- X_model_matrix_1_sub[, !colnames(X_model_matrix_1_sub) %in% correlated]
  }


  YXZ_2_sub <- data.frame(target = ph_sub)
  xnames_1 <- colnames(X_model_matrix_1_sub)
  YXZ_2_sub[, xnames_1] <- X_model_matrix_1_sub
  YXZ_2_sub[, znames_1] <- Z2[!missind, , drop = FALSE]
  YXZ_2_sub[, "clID"] <- clID[!missind]

  tmp_2_all <- tmp_0_all[, colnames(YXZ_2_sub), drop = FALSE]

  # -------------- calling the gibbs sampler to get imputation parameters----
  fixformula <- stats::formula(paste("target~ - 1 + ",
                                       paste(xnames_1, sep = "", collapse = "+")))

  randformula <- stats::formula(paste("~us( - 1 + ", paste(znames_1, sep = "", collapse = "+"),
                                 "):clID", sep = ""))

  number_fix_parameters <- ncol(X_model_matrix_1_sub)

  # Get the number of random effects variables
  number_random_effects <- length(znames_1)
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
                                       data = YXZ_2_sub,
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

  fix.eff.imp <- matrix(xdraws[select.record, ], nrow = ncol(X_model_matrix_1_sub))

  sigma.y.imp <- sqrt(variancedraws[select.record, ncol(variancedraws)])

  latent <- stats::rnorm(n, as.matrix(tmp_2_all[, xnames_1, drop = FALSE]) %*% fix.eff.imp +
                      apply(Z2 * rand.eff.imp[clID,], 1, sum), sigma.y.imp)

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
