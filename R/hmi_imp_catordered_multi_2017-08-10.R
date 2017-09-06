# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'
#
# MS: All other shortcuts are shown when pressing 'Alt + Shift + K'
# MS: Hinweis: es wird zu jeder R-Datei im Projekordner\R, die eine Documentation hat,
#Mit devtools::document() auch eine .RD Datei erstellt.


#' The function for hierarchical imputation of categorical variables.
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
#' @return A n x 1 data.frame with the original and imputed values as a factor.
imp_orderedcat_multi <- function(y_imp,
                      X_imp,
                      Z_imp,
                      clID,
                      nitt = 3000,
                      thin = 10,
                      burnin = 1000){

  if(!is.factor(y_imp)){
    warning("We suggest to make all your categorical variables to be a factor.")
    y_imp <- as.factor(y_imp)
  }

  # -----------------------------preparing the data ------------------
  # -- standardise the covariates in X (which are numeric and no intercept)
  # remove excessive variables
  X_imp <- remove_excessives(X_imp)

  X_imp_stand <- stand(X_imp)

  # -- standardise the covariates in Z (which are numeric and no intercept)
  Z_imp_stand <- stand(Z_imp)

  n <- nrow(X_imp_stand)
  missind <- is.na(y_imp)

  # Get the number of random effects variables
  n.par.rand <- ncol(Z_imp_stand)
  length.alpha <- length(table(clID)) * n.par.rand

  number_clusters <- length(table(clID))

  #starting model
  ph <- sample_imp(y_imp)[, 1]
  tmp_0_all <- data.frame(target = ph)
  xnames_0 <- paste("X", 1:ncol(X_imp_stand), sep = "")
  tmp_0_all[xnames_0] <- X_imp_stand
  tmp_0_sub <- tmp_0_all[!missind, , drop = FALSE]


  types <- array(dim = ncol(tmp_0_all))
  for(i in 1:length(types)) types[i] <- get_type(tmp_0_all[, i])
  #polr needs an intercept variable. So we have to include one to the model formuly,
  #but not if the data already have an intercept variable.

  tmp_0_all_nointercept <- tmp_0_all[, types != "intercept", drop = FALSE]
  tmp_0_sub_nointercept <- tmp_0_sub[, types != "intercept", drop = FALSE]

  reg_1_all <- MASS::polr(target ~ 1 + ., data =  tmp_0_all_nointercept, method = "probit")
  reg_1_sub <- MASS::polr(target ~ 1 + ., data =  tmp_0_sub_nointercept, method = "probit")



  X_model_matrix_1_all <- stats::model.matrix(reg_1_all)
  xnames_1 <- paste("X", 1:ncol(X_model_matrix_1_all), sep = "")
  znames <- paste("Z", 1:ncol(Z_imp_stand), sep = "")

  #remove unneeded variables
  unneeded <- is.na(stats::coefficients(reg_1_sub))
  xnames_2 <- xnames_1[!unneeded]

  tmp_2_all <- data.frame(target = as.factor(y_imp))
  X_model_matrix_2_all <- X_model_matrix_1_all[, !unneeded, drop = FALSE]
  tmp_2_all[, xnames_2] <- X_model_matrix_2_all
  tmp_2_all[, znames] <- Z_imp_stand
  tmp_2_all[, "ClID"] <- clID
  tmp_2_sub <- tmp_2_all[!missind, , drop = FALSE]

  # -------------- calling the gibbs sampler to get imputation parameters----

  fixformula <- stats::formula(paste("target~ - 1 + ",
                                paste(xnames_2, sep = "", collapse = "+")))

  randformula <- stats::formula(paste("~us( - 1 + ", paste(znames, sep = "", collapse = "+"),
                                 "):ClID", sep = ""))



  number_fix_parameters <- ncol(X_model_matrix_2_all)

  # Get the number of random effects variables
  number_random_effects <- length(znames)
  number_random_parameters <- number_random_effects
  #Fix residual variance R at 1
  # cf. http://stats.stackexchange.com/questions/32994/what-are-r-structure-g-structure-in-a-glmm

  J <- length(table(y_imp)) #number of categories
  #priors from https://stat.ethz.ch/pipermail/r-sig-mixed-models/2012q2/018194.html

  prior <- list(R = list(V = 1, nu = 0.002),
       G = list(G1 = list(V = diag(number_random_effects), nu = 0.002)))

  MCMCglmm_draws <- MCMCglmm::MCMCglmm(fixed = fixformula,
                                       random = randformula,
                                       data = tmp_2_sub,
                                       family = "ordinal",
                                       verbose = FALSE, pr = TRUE,
                                       prior = prior,
                                       saveX = TRUE, saveZ = TRUE,
                                       nitt = 1e5,
                                       thin = 1e3,
                                       burnin = 1e4)

  pointdraws <- MCMCglmm_draws$Sol

  xdraws <- pointdraws[, 1:number_fix_parameters, drop = FALSE]
  zdraws <- pointdraws[, number_fix_parameters +
                         1:(number_random_parameters * number_clusters), drop = FALSE]

  variancedraws <- MCMCglmm_draws$VCV

  number_of_draws <- nrow(pointdraws)
  select.record <- sample(1:number_of_draws, 1, replace = TRUE)

  # -------------------- drawing samples with the parameters from the gibbs sampler --------

  #now decide in with category a person falls acording to threshold.



  ###start imputation



  rand.eff.imp <- matrix(zdraws[select.record,],
                           ncol = number_random_effects) # cluster specific effects

  fix.eff.imp <- matrix(xdraws[select.record, ], nrow = ncol(X_model_matrix_2_all))

  sigma.y.imp <- sqrt(variancedraws[select.record, ncol(variancedraws)])

  latent <- stats::rnorm(n, X_model_matrix_2_all %*% fix.eff.imp +
                      apply(Z_imp_stand * rand.eff.imp[clID,], 1, sum), sigma.y.imp)

    #cutpoint: the first cutpoint is set to 0 by MCMCglmm
    #(see https://stat.ethz.ch/pipermail/r-sig-mixed-models/2013q1/020030.html)

  cps <- c(min(latent), 0, MCMCglmm_draws$CP[select.record,], max(latent))

  y_temp <- factor(levels(y_imp)[as.numeric(cut(latent, breaks = cps, include.lowest = TRUE))],
                ordered = TRUE)


  y_ret <- factor(ifelse(missind, as.character(y_temp),
                                   as.character(y_imp)),
                            ordered = TRUE)

  # --------- returning the imputed data --------------
  return(data.frame(y_imp = y_ret))
}
