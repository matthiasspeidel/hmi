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
#' @return A list with 1. 'y_ret' the n x 1 data.frame with the original and imputed values.
#' 2. 'Sol' the Gibbs-samples for the fixed effects parameters.
#' 3. 'VCV' the Gibbs-samples for variance parameters.
imp_binary_multi <- function(y_imp,
                      X_imp,
                      Z_imp,
                      clID,
                      nitt = 22000,
                      burnin = 2000,
                      thin = 20){

  # ----------------------------- preparing the y data ------------------
  # stransform y_imp into a real binary with only zeros and ones (and NAs).
  first_possibility <- utils::head(sort(y_imp), n = 1)
  second_possibility <- utils::tail(sort(y_imp), n = 1)
  y_binary <- data.frame(y = factor(y_imp, labels = c(0, 1)))


  # If one category has less then two observations, no binary model can be estimated.
  # So the imputation routines has to stop.
  if(min(table(y_binary)) < 2){
    stop("A binary (or maybe a semicontinuous) variable has less than two observations in one category.
         Consider removing this variable
         (or in the case of a semicontinuous variable, to specify it as continuous in the list_of_types (see ?hmi)).")
  }

  # ----------------------------- preparing the X and Z data ------------------

  # remove excessive variables
  X_imp <- cleanup(X_imp)

  # standardise the covariates in X (which are numeric and no intercept)
  X_imp_stand <- stand(X_imp)

  # -- standardise the covariates in Z (which are numeric and no intercept)
  Z_imp <- cleanup(Z_imp)
  Z_imp_stand <- stand(Z_imp)


  missind <- is.na(y_binary)
  n <- nrow(X_imp_stand)

  # Get the number of random effects variables
  n.par.rand <- ncol(Z_imp_stand)
  length.alpha <- length(table(clID)) * n.par.rand

  # We need two design matrices. One for the model to get the imputation parameters
  # which is based on the observed values only (obs).
  # And secondly a design matrix for E(y_mis|X_mis) based on all observations (all).
  # The later is more or less only the technical framework to get the imputed values.

  #define a place holder (ph)
  ph <- sample_imp(y_binary[, 1])[, 1]


  tmp_0_sub <- data.frame(target = ph, X_imp_stand)[!missind, , drop = FALSE]
  tmp_0_all <- data.frame(target = ph, X_imp_stand)


  X_model_matrix_1_sub <- stats::model.matrix(target ~ 0 + ., data = tmp_0_sub)
  X_model_matrix_1_all <- stats::model.matrix(target ~ 0 + ., data = tmp_0_all)
  colnames(X_model_matrix_1_sub) <- gsub("`", "", colnames(X_model_matrix_1_sub))
  colnames(X_model_matrix_1_all) <- gsub("`", "", colnames(X_model_matrix_1_all))


  # remove unneeded variables/categories from X_model_matrix_1
  # model to determine unnneeded variables
  reg_1_sub <- stats::glm(target ~ 0 + ., data = tmp_0_sub,
                          family = stats::binomial(link = "logit"))

  #data, where the needed variables should be stored
  tmp_2_all <- data.frame(target = ph)

  xnames_1 <- paste("X", 1:ncol(X_model_matrix_1_all), sep = "")
  znames_1 <- paste("Z", 1:ncol(Z_imp_stand), sep = "")

  xnames_2 <- xnames_1[!is.na(stats::coefficients(reg_1_sub))]
  znames_2 <- znames_1

  tmp_2_all[, xnames_2] <- X_model_matrix_1_all[, !is.na(stats::coefficients(reg_1_sub)),
                                                drop = FALSE]
  tmp_2_all[, znames_2] <- Z_imp_stand
  tmp_2_all[, "ClID"] <- clID

  tmp_2_sub <- tmp_2_all[!missind, , drop = FALSE]


  glmfixformula_2 <- stats::formula(paste("target ~ 0 +",
                                         paste(xnames_2, collapse = "+"), sep = ""))

  reg_2_all <- stats::glm(glmfixformula_2, data = tmp_2_all,
                          family = stats::binomial(link = "logit"))


  X_model_matrix_2_all <- stats::model.matrix(reg_2_all)

  # -------------- calling the gibbs sampler to get imputation parameters ----

  fixformula <- stats::formula(paste("target~", paste(xnames_2, collapse = "+"), "- 1",
                                       sep = ""))
  randformula <- stats::as.formula(paste("~us(0+", paste(znames_2, collapse = "+"), "):ClID",
                                           sep = ""))

  #Fix residual variance R at 1
  # cf. http://stats.stackexchange.com/questions/32994/what-are-r-structure-g-structure-in-a-glmm
  prior <- list(R = list(V = 1, fix = TRUE),
                G = list(G1 = list(V = diag(n.par.rand), nu = 0.002)))

  MCMCglmm_draws <- MCMCglmm::MCMCglmm(fixed = fixformula,
                                       random = randformula,
                                       data = tmp_2_sub,
                                       family = "categorical",
                                       verbose = FALSE, pr = TRUE, prior = prior,
                                       saveX = TRUE, saveZ = TRUE,
                                       nitt = nitt,
                                       thin = thin,
                                       burnin = burnin)

  # correction. see:
  # http://stats.stackexchange.com/questions/32994/what-are-r-structure-g-structure-in-a-glmm
  k <- ((16*sqrt(3))/(15*pi))^2

  pointdraws <- MCMCglmm_draws$Sol /sqrt(1 + k)
  xdraws <- pointdraws[, 1:ncol(X_model_matrix_2_all), drop = FALSE]
  zdraws <- pointdraws[, ncol(X_model_matrix_2_all) + 1:length.alpha, drop = FALSE]
  variancedraws <- MCMCglmm_draws$VCV

  number_of_draws <- nrow(pointdraws)
  select_record <- sample(1:number_of_draws, size = 1)

  # -------------------- drawing samples with the parameters from the gibbs sampler --------


  linkfunction <- boot::inv.logit


  ###start imputation

  rand_eff_imp <- matrix(zdraws[select_record,],
                           ncol = n.par.rand)

  fix_eff_imp <- matrix(xdraws[select_record, ], nrow = ncol(X_model_matrix_2_all))

  sigma_y_imp <- sqrt(variancedraws[select_record, ncol(variancedraws)])

  linearpredictor <- stats::rnorm(n, X_model_matrix_2_all %*% fix_eff_imp +
                      apply(Z_imp_stand * rand_eff_imp[clID,], 1, sum), 0*sigma_y_imp)

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
