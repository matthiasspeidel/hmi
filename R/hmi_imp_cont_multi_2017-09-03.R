#' The function for hierarchical imputation of continuous variables.
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
#' @return A n x 1 matrix with the original and imputed values.
imp_cont_multi <- function(y_imp,
                      X_imp,
                      Z_imp,
                      clID,
                      nitt = 3000,
                      thin = 10,
                      burnin = 1000){

  # -----------------------------preparing the data ------------------
  # -- standardise the covariates in X (which are numeric and no intercept)
  # ----------------------------- preparing the X and Z data ------------------

  # remove excessive variables
  X_imp <- remove_excessives(X_imp)

  # standardise the covariates in X (which are numeric and no intercept)
  X_imp_stand <- stand(X_imp)

  # -- standardise the covariates in Z (which are numeric and no intercept)

  Z_imp_stand <- stand(Z_imp)


  #define a place holder (ph)
  ph <- sample_imp(y_imp)[, 1]

  y_mean <- mean(ph, na.rm = TRUE)
  y_sd <- stats::sd(ph, na.rm = TRUE)

  ph_stand <- (ph - y_mean)/y_sd + 1
  #Z may contain factor variables. Later we need Z to contain numeric variables,
  # so we have to change those variables.
  tmp <- data.frame(target = ph_stand, Z_imp_stand)
  #remove intercept variable
  tmp <- tmp[, get_type(tmp) != "intercept", drop = FALSE]

  Z_imp_stand2 <- stats::model.matrix(stats::lm("target ~ 1 + .", data = tmp))

  missind <- is.na(y_imp)
  n <- nrow(X_imp_stand)

  tmp_0_sub <- data.frame(target = ph_stand, X_imp_stand)[!missind, , drop = FALSE]
  tmp_0_all <- data.frame(target = ph_stand, X_imp_stand)


  X_model_matrix_1_sub <- stats::model.matrix(target ~ 0 + ., data = tmp_0_sub)
  X_model_matrix_1_all <- stats::model.matrix(target ~ 0 + ., data = tmp_0_all)
  colnames(X_model_matrix_1_sub) <- gsub("`", "", colnames(X_model_matrix_1_sub))
  colnames(X_model_matrix_1_all) <- gsub("`", "", colnames(X_model_matrix_1_all))

  # remove unneeded variables/categories from X_model_matrix_1
  # model to determine unnneeded variables
  reg_1_sub <- stats::lm(target ~ 0 + ., data = tmp_0_sub)
  unneeded <- is.na(stats::coefficients(reg_1_sub))

  #data, where the needed variables should be stored


  xnames_1 <- paste("X", 1:ncol(X_model_matrix_1_all), sep = "")
  znames_1 <- paste("Z", 1:ncol(Z_imp_stand2), sep = "")

  xnames_2 <- xnames_1[!unneeded]
  znames_2 <- znames_1
  X_model_matrix_2_all <- X_model_matrix_1_all[, !unneeded, drop = FALSE]
  tmp_2_all <- data.frame(target = ph_stand)
  tmp_2_all[, xnames_2] <- X_model_matrix_2_all
  tmp_2_all[, znames_2] <- Z_imp_stand2
  tmp_2_all[, "ClID"] <- clID

  tmp_2_sub <- tmp_2_all[!missind, , drop = FALSE]
  # -------------- calling the gibbs sampler to get imputation parameters----

  tmp_3_sub <- tmp_2_sub
  intercept_variables <- get_type(tmp_3_sub) == "intercept"
  tmp_3_sub <- tmp_2_sub[, !intercept_variables, drop = FALSE]
  xnames_3 <- xnames_2[!xnames_2 %in% names(intercept_variables)[intercept_variables]]
  znames_3 <- znames_2[!znames_2 %in% names(intercept_variables)[intercept_variables]]
  fixformula <- stats::formula(paste("target~ 1+ ", paste(xnames_3, collapse = "+"), sep = ""))
  randformula <- stats::as.formula(paste("~us(1+", paste(znames_3, collapse = "+"), "):ClID", sep = ""))


  prior <- list(R = list(V = 1, nu = 0.002), # alternatice: R = list(V = 1e-07, nu = -2)
                G = list(G1 = list(V = diag(ncol(Z_imp_stand2)), nu = 0.002)))


  MCMCglmm_draws <- MCMCglmm::MCMCglmm(fixformula, random = randformula, data = tmp_3_sub,
                       verbose = FALSE, pr = TRUE, prior = prior,
                       saveX = TRUE, saveZ = TRUE,
                       nitt = 10000,
                       thin = 100,
                       burnin = 5000)

  # Get the number of random effects variables
  n.par.rand <- ncol(Z_imp_stand2)
  ncluster <- length(table(tmp_2_sub$ClID))
  length.alpha <- ncluster * n.par.rand

  pointdraws <- MCMCglmm_draws$Sol
  xdraws <- pointdraws[, 1:ncol(X_model_matrix_2_all), drop = FALSE]
  #If a cluster cannot has random effects estimates, because there too few observations,
  #we make them 0.
  empty_cluster <- which(table(tmp_2_sub$ClID) == 0)
  zdraws_pre <- pointdraws[, (ncol(X_model_matrix_2_all) + 1):ncol(pointdraws), drop = FALSE]
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
  select.record <- sample(1:number_of_draws, 1, replace = TRUE)

  # -------------------- drawing samples with the parameters from the gibbs sampler --------
  ###start imputation


  rand.eff.imp <- matrix(zdraws[select.record, ], ncol = n.par.rand)

  fix.eff.imp <- matrix(xdraws[select.record, ], nrow = ncol(X_model_matrix_2_all))

  sigma.y.imp <- sqrt(variancedraws[select.record, ncol(variancedraws)])

  y_temp <- stats::rnorm(n, X_model_matrix_2_all %*% fix.eff.imp +
                      apply(Z_imp_stand2 * rand.eff.imp[clID,], 1, sum), sd = sigma.y.imp)

  y_ret <- matrix(ifelse(missind, (y_temp - 1) * y_sd + y_mean, y_imp), ncol = 1)


  # --------- returning the imputed data --------------
  return(y_ret)
}


# Generate documentation with devtools::document()
# Build package with devtools::build() and devtools::build(binary = TRUE) for zips
