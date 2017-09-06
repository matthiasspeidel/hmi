#' The function for hierarchical imputation of variables with count data.
#'
#' The function is called by the wrapper.
#' @param y_imp A Vector with the variable to impute.
#' @param X_imp A data.frame with the fixed effects variables.
#' @param Z_imp A data.frame with the random effects variables.
#' @param clID A vector with the cluster ID.
#' @param nitt An integer defining number of MCMC iterations (see MCMCglmm).
#' @param thin An integer defining the thinning interval (see MCMCglmm).
#' @param burnin An integer defining the percentage of draws from the gibbs sampler
#' that should be discarded as burn in (see MCMCglmm).
#' @return A n x 1 data.frame with the original and imputed values.
imp_count_multi <- function(y_imp,
                      X_imp,
                      Z_imp,
                      clID,
                      nitt = 3000,
                      thin = 10,
                      burnin = 1000){
  # ----------------------------- preparing the X and Z data ------------------

  # remove excessive variables
  X_imp <- remove_excessives(X_imp)

  # standardise the covariates in X (which are numeric and no intercept)
  X_imp_stand <- stand(X_imp)

  # -- standardise the covariates in Z (which are numeric and no intercept)

  Z_imp_stand <- stand(Z_imp)

  #the missing indactor indicates, which values of y are missing.
  missind <- is.na(y_imp)
  n <- nrow(X_imp_stand)
  number_clusters <- length(table(clID))

  # Get the number of random effects variables
  n.par.rand <- ncol(Z_imp_stand)
  length.alpha <- length(table(clID)) * n.par.rand

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
  znames <- paste("Z", 1:ncol(Z_imp_stand), sep = "")

  #remove unneeded variables
  unneeded <- is.na(stats::coefficients(reg_1_sub))
  xnames_2 <- xnames_1[!unneeded]

  tmp_2_all <- data.frame(target = y_imp)
  X_model_matrix_2_all <- X_model_matrix_1_all[, !unneeded, drop = FALSE]

  tmp_2_all[, xnames_2] <- X_model_matrix_2_all
  # note: the [] part is only relevant if ncol(X_model_matrix) == 1
  # (see http://stackoverflow.com/questions/40872034/what-happens-when-a-data-frame-gets-new-columns)
  tmp_2_all[, znames] <- Z_imp_stand[, 1:ncol(Z_imp_stand)]
  tmp_2_all[, "ClID"] <- clID

  tmp_2_sub <- tmp_2_all[!missind, , drop = FALSE]

  # -------------- calling the gibbs sampler to get imputation parameters----

  fixformula <- stats::formula(paste("target~", paste(xnames_2, collapse = "+"), "- 1",
                                     sep = ""))
  randformula <- stats::as.formula(paste("~us(", paste(znames, collapse = "+"), "):ClID",
                                         sep = ""))


  prior <- list(R = list(V = 1e-07, nu = -2),
                G = list(G1 = list(V = diag(ncol(Z_imp_stand)), nu = 0.002)))


  MCMCglmm_draws <- MCMCglmm::MCMCglmm(fixformula, random = randformula,
                                       data = tmp_2_sub,
                                       verbose = FALSE, pr = TRUE, prior = prior,
                                       family = "poisson",
                                       saveX = TRUE, saveZ = TRUE,
                                       nitt = nitt,
                                       thin = thin,
                                       burnin = burnin)


  pointdraws <- MCMCglmm_draws$Sol
  xdraws <- pointdraws[, 1:ncol(X_model_matrix_2_all), drop = FALSE]
  zdraws <- pointdraws[, ncol(X_model_matrix_2_all) + 1:length.alpha, drop = FALSE]
  variancedraws <- MCMCglmm_draws$VCV
  # the last column contains the variance (not standard deviation) of the residuals

  number_of_draws <- nrow(pointdraws)
  select.record <- sample(1:number_of_draws, 1, replace = TRUE)

  # -------------------- drawing samples with the parameters from the gibbs sampler --------
  ###start imputation


  rand.eff.imp <- matrix(zdraws[select.record,],
                           ncol = n.par.rand)

  fix.eff.imp <- matrix(xdraws[select.record, ], nrow = ncol(X_model_matrix_2_all))

  sigma.y.imp <- sqrt(variancedraws[select.record, ncol(variancedraws)])

  lambda <- exp(stats::rnorm(n, X_model_matrix_2_all %*% fix.eff.imp +
                      apply(Z_imp_stand * rand.eff.imp[clID,], 1, sum), sigma.y.imp))

  y_proposal <- stats::rpois(n, lambda)
  y_tmp <- ifelse(missind, y_proposal, y_imp)


  # --------- returning the imputed data --------------
  return(data.frame(y_imp = y_tmp))

}


# Generate documentation with devtools::document()
# Build package with devtools::build() and devtools::build(binary = TRUE) for zips
