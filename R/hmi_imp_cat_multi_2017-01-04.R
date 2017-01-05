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
#' The function is called by the wrapper and relies on \code{MCMCglmm}.\cr
#' While in the single level function (\code{imp_cat_single}) we used regression trees
#' to impute data, here we run a multilevel multinomial model.
#' The basic idea is that for each category of the target variable (expect the reference category)
#' a own formula is set up, saying for example that the chances to end up in category
#' j increase with increasing X5. So there is a own regression coefficient beta_5_j present.
#' In a multilevel setting, this regression coefficient beta_5_j might be different for
#' different clusters so for cluster 1 it would be beta_5_j_1 = beta_5_j + u_5_1 and for
#' cluster 27 beta_5_j_27 = beta_5_j + u_5_27. This also leads to own random effect covariance
#' matrices for each category. Or, if you want to have all random effect variance parameters
#' in one matrix: a (very large) matrix where not for example only the random intercepts variance
#' and random slopes variance and their covariance is present. Instead there is even a
#' covariance between the random slopes in category 2 and the random intercepts in category 4.
#' For simplicity these covariances are set to be 0.
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
imp_cat_multi <- function(y_imp_multi,
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

  n <- nrow(X_imp_multi_stand)
  tmp0 <- data.frame(y_tmp = rnorm(n), X_imp_multi_stand, cl.id = clID)
  tmp0_formula <- formula(paste("y_tmp ~ 0 + ", paste(colnames(X_imp_multi_stand), collapse = "+")))

  X_model_matrix <- model.matrix(tmp0_formula, data = tmp0)

  # Remove ` from the variable names
  colnames(X_model_matrix) <- gsub("`", "", colnames(X_model_matrix))



  # -- standardise the covariates in Z (which are numeric and no intercept)
  need_stand <- apply(Z_imp_multi, 2, get_type) == "cont"
  Z_imp_multi_stand <- Z_imp_multi

  number_clusters <- length(table(clID))


  # -------------- calling the gibbs sampler to get imputation parameters----

  tmp <- data.frame(target = y_imp_multi)
  xnames <- paste("X", 1:ncol(X_model_matrix), sep = "")
  znames <- paste("Z", 1:ncol(Z_imp_multi_stand), sep = "")
  tmp[, xnames] <- X_model_matrix[, 1:ncol(X_model_matrix)]
  # note: the [] part is only relevant if ncol(X_model_matrix) == 1
  # (see http://stackoverflow.com/questions/40872034/what-happens-when-a-data-frame-gets-new-columns)
  tmp[, znames] <- Z_imp_multi_stand[, 1:ncol(Z_imp_multi_stand)]
  tmp[, "ClID"] <- clID

  fixformula <- formula(paste("target~ - 1 + ",
                                paste(xnames, ":trait", sep = "", collapse = "+"), sep = ""))

  randformula <- formula(paste("~us( - 1 + ", paste(znames, ":trait", sep = "", collapse = "+"),
                                 "):ClID", sep = ""))

  J <- length(table(y_imp_multi)) #number of categories
  number_fix_parameters <- ncol(X_model_matrix) * (J-1)

  # Get the number of random effects variables
  number_random_effects <- length(znames)
  number_random_parameters <- number_random_effects * (J - 1)
  #Fix residual variance R at 1
  # cf. http://stats.stackexchange.com/questions/32994/what-are-r-structure-g-structure-in-a-glmm


  J_matrix <- array(1, dim = c(J, J) -1) # matrix of ones
  I_matrix <- diag(J - 1) #identiy matrix

  IJ <- (I_matrix + J_matrix)/J # see Hadfields Course notes p 97
  prior <- list(R = list(V = IJ, fix = TRUE),
                G = list(G1 = list(V = diag(number_random_parameters), nu = 2)))

  MCMCglmm_draws <- MCMCglmm::MCMCglmm(fixed = fixformula,
                                       random = randformula,
                                       rcov = ~us(trait):units,
                                       data = tmp,
                                       family = "categorical",
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


  number_of_draws <- nrow(pointdraws)
  select.record <- sample(1:number_of_draws, M, replace = TRUE)

  # -------------------- drawing samples with the parameters from the gibbs sampler --------
  #now generate new P(Y = A|x * beta) = x*beta/(1+ sum(exp(x*beta))) etc.

  #set up random intercepts and slopes


  exp_beta <- array(dim = c(n, J - 1))

  fix_eff_sample <- xdraws[select.record[1], ]
  rand_eff_sample <- zdraws[select.record[1], ]
  for(k in 1:ncol(exp_beta)){

    rand_betas <- matrix(rand_eff_sample[(k - 1) * number_clusters +
                     c(1:number_clusters, 1:number_clusters  + (J - 1) * number_clusters)], nrow = number_clusters)

    exp_beta[, k] <-
      as.matrix(exp(as.matrix(tmp[, xnames, drop = FALSE]) %*%
                      fix_eff_sample[seq(0, J*ncol(exp_beta), by = J - 1)[1:ncol(exp_beta)] + k] +#fix effects
      rowSums(tmp[, znames, drop = FALSE] * rand_betas[clID, , drop = FALSE])))#random effects

    #explanation for the fixed effects part:
    #MCMCglmm_draws$Sol is ordered in the following way: beta_1 for category 1, beta_1 for category_2
    #... beta 1 for category_k, beta_2 for category_1, beta_2 for category_2 ...
    #so we have to skip some values as we proceed by category and not beta.
  }

  y_imp <- data.frame(matrix(nrow = n, ncol = M))
  ###start imputation
  for (j in 1:M){

    y_temp <- array(dim = n)


    for(i in 1:n){

      # enshure that the reference category is in the correct position
      mytmp_prob <- exp_beta[i, ]/(1 + sum(exp_beta[i, ]))
      my_prob <- c(1 - sum(mytmp_prob), mytmp_prob) #the first category is the reference category in MCMCglmm

      y_temp[i] <- sample(levels(y_imp_multi), size = 1, replace = TRUE, prob = my_prob)

    }

    y_imp[, j] <- as.factor(ifelse(is.na(y_imp_multi), y_temp, as.character(y_imp_multi)))
  }
  colnames(y_imp) <- paste("imp", 1:M, sep = "_")

  # --------- returning the imputed data --------------
  return(y_imp)
}


# Generate documentation with devtools::document()
# Build package with devtools::build() and devtools::build(binary = TRUE) for zips
