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
#' @param y_imp_multi A Vector with the variable to impute.
#' @param X_imp_multi A data.frame with the fixed effects variables.
#' @param Z_imp_multi A data.frame with the random effects variables.
#' @param clID A vector with the cluster ID.
#' @param model_formula A \code{\link[stats]{formula}} used for the analysis model.
#' @param M An integer defining the number of imputations that should be made.
#' @param nitt An integer defining number of MCMC iterations (see MCMCglmm).
#' @param thin An integer defining the thinning interval (see MCMCglmm).
#' @param burnin An integer defining the percentage of draws from the gibbs sampler
#' that should be discarded as burn in (see MCMCglmm).
#' @return A n x M matrix. Each column is one of M imputed y-variables.
imp_orderedcat_multi <- function(y_imp_multi,
                      X_imp_multi,
                      Z_imp_multi,
                      clID,
                      model_formula,
                      M = 10,
                      nitt = 3000,
                      thin = 10,
                      burnin = 1000){

  if(!is.factor(y_imp_multi)){
    warning("We suggest to make all your categorical variables to be a factor.")
    y_imp_multi <- as.factor(y_imp_multi)
  }

  # -----------------------------preparing the data ------------------
  # -- standardise the covariates in X (which are numeric and no intercept)
  need_stand <- apply(X_imp_multi, 2, get_type) == "cont"
  X_imp_multi_stand <- X_imp_multi
  X_imp_multi_stand[, need_stand] <- scale(X_imp_multi[, need_stand])

  #generate model.matrix (from the class matrix)
  n <- nrow(X_imp_multi_stand)
  X_model_matrix <- stats::model.matrix(stats::rnorm(n) ~ 0 + ., data = X_imp_multi_stand)
  # Remove ` from the variable names
  colnames(X_model_matrix) <- gsub("`", "", colnames(X_model_matrix))

  # -- standardise the covariates in Z (which are numeric and no intercept)
  need_stand <- apply(Z_imp_multi, 2, get_type) == "cont"
  Z_imp_multi_stand <- Z_imp_multi
  Z_imp_multi_stand[, need_stand] <- scale(Z_imp_multi[, need_stand])

  # Get the number of random effects variables
  n.par.rand <- ncol(Z_imp_multi_stand)
  length.alpha <- length(table(clID)) * n.par.rand


  # -------------- calling the gibbs sampler to get imputation parameters----


  n <- length(y_imp_multi)
  lmstart <- stats::lm(stats::rnorm(n) ~ 0 +., data = X_imp_multi)

  X_model_matrix_1 <- stats::model.matrix(lmstart)
  xnames_1 <- paste("X", 1:ncol(X_model_matrix_1), sep = "")
  znames <- paste("Z", 1:ncol(Z_imp_multi_stand), sep = "")

  tmp_1 <- data.frame(y = stats::rnorm(n))
  tmp_1[, xnames_1] <- X_model_matrix_1

  reg_1 <- stats::lm(y ~ 0 + . , data = tmp_1)

  blob <- y_imp_multi
  tmp_2 <- data.frame(target = blob)

  xnames_2 <- xnames_1[!is.na(stats::coefficients(reg_1))]
  X_model_matrix_2 <- X_model_matrix_1[, !is.na(stats::coefficients(reg_1)), drop = FALSE]
  tmp_2[, xnames_2] <- X_model_matrix_2
  tmp_2[, znames] <- Z_imp_multi_stand
  tmp_2[, "ClID"] <- clID


  number_clusters <- length(table(clID))


  # -------------- calling the gibbs sampler to get imputation parameters----


  fixformula <- stats::formula(paste("target~ - 1 + ",
                                paste(xnames_2, sep = "", collapse = "+")))

  randformula <- stats::formula(paste("~us( - 1 + ", paste(znames, sep = "", collapse = "+"),
                                 "):ClID", sep = ""))


  number_fix_parameters <- ncol(X_model_matrix_2)

  # Get the number of random effects variables
  number_random_effects <- length(znames)
  number_random_parameters <- number_random_effects
  #Fix residual variance R at 1
  # cf. http://stats.stackexchange.com/questions/32994/what-are-r-structure-g-structure-in-a-glmm

  J <- length(table(y_imp_multi)) #number of categories
  #priors from https://stat.ethz.ch/pipermail/r-sig-mixed-models/2012q2/018194.html

  prior <- list(R = list(V = 1, fix = TRUE),
       G = list(G1 = list(V = diag(number_random_effects), nu = 2)))

  MCMCglmm_draws <- MCMCglmm::MCMCglmm(fixed = fixformula,
                                       random = randformula,
                                       data = tmp_2,
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
  select.record <- sample(1:number_of_draws, M, replace = TRUE)

  # -------------------- drawing samples with the parameters from the gibbs sampler --------

  #now decide in with category a person falls acording to threshold.


  y_imp <- data.frame(matrix(nrow = n, ncol = M))
  ###start imputation
  for (j in 1:M){


    rand.eff.imp <- matrix(zdraws[select.record[j],],
                           ncol = number_random_effects) # cluster specific effects

    fix.eff.imp <- matrix(xdraws[select.record[j], ], nrow = ncol(X_model_matrix_2))

    sigma.y.imp <- sqrt(variancedraws[select.record[j], ncol(variancedraws)])

    latent <- stats::rnorm(n, X_model_matrix_2 %*% fix.eff.imp +
                      apply(Z_imp_multi_stand * rand.eff.imp[clID,], 1, sum), sigma.y.imp)

    #cutpoint: the first cutpoint is set to 0 by MCMCglmm
    #(see https://stat.ethz.ch/pipermail/r-sig-mixed-models/2013q1/020030.html)

    cps <- c(min(latent), 0, MCMCglmm_draws$CP[select.record[j],], max(latent))

    y_temp <- factor(levels(y_imp_multi)[as.numeric(cut(latent, breaks = cps, include.lowest = TRUE))],
                ordered = TRUE)


    y_imp[, j] <- factor(ifelse(is.na(y_imp_multi),
                                   as.character(y_temp),
                                   as.character(y_imp_multi)),
                            ordered = TRUE)
  }
  colnames(y_imp) <- paste("imp", 1:M, sep = "_")

  # --------- returning the imputed data --------------
  return(y_imp)
}


# Generate documentation with devtools::document()
# Build package with devtools::build() and devtools::build(binary = TRUE) for zips
