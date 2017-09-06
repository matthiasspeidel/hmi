#' The function for hierarchical imputation of categorical variables.
#'
#' The function is called by the wrapper and relies on \code{MCMCglmm}.\cr
#' While in the single level function (\code{imp_cat_single}) we used regression trees
#' to impute data, here we run a multilevel multinomial model.
#' The basic idea is that for each category of the target variable (expect the reference category)
#' an own formula is set up, saying for example that the chances to end up in category
#' j increase with increasing X5. So there is an own regression coefficient \eqn{beta_{5,j}} present.
#' In a multilevel setting, this regression coefficient \eqn{beta_{5,j}} might be different for
#' different clusters: for cluster 27 it would be \eqn{beta_{5,j,27} = beta_{5,j} + u_{5,27}}.
#' This also leads to own random effect covariance matrices for each category.
#' All those random effect variance parameters can be collected
#' in one (quite large) covariance matrix where (for example)
#' not only the random intercepts variance and random slopes variance and their covariance
#' is present. Instead, there is even a covariance between the random slopes in category s
#' and the random intercepts in category p. Beside the difficulties in interpretation,
#' these covariances have shown to be numerically instable so they are set to be 0.
#' @param y_imp A Vector with the variable to impute.
#' @param X_imp A data.frame with the fixed effects variables.
#' @param Z_imp A data.frame with the random effects variables.
#' @param clID A vector with the cluster ID.
#' @param nitt An integer defining number of MCMC iterations (see MCMCglmm).
#' @param thin An integer defining the thinning interval (see MCMCglmm).
#' @param burnin An integer defining the percentage of draws from the gibbs sampler
#' that should be discarded as burn in (see MCMCglmm).
#' @return A n x 1 data.frame with the original and imputed values.
imp_cat_multi <- function(y_imp,
                      X_imp,
                      Z_imp,
                      clID,
                      nitt = 1e5,
                      thin = 1e3,
                      burnin = 1e4){

  if(min(table(y_imp)) < 2) {
    stop("Too few observations per category in a categorical target variable.")
  }

  if(!is.factor(y_imp)){
    warning("We suggest to make all your categorical variables to be a factor.")
    y_imp <- as.factor(y_imp)
  }

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


  #starting model
  ph <- sample_imp(y_imp)[, 1]
  tmp_0_all <- data.frame(target = ph)
  xnames_0 <- paste("X", 1:ncol(X_imp_stand), sep = "")
  tmp_0_all[xnames_0] <- X_imp_stand
  tmp_0_sub <- tmp_0_all[!missind, , drop = FALSE]

  reg_1_all <- nnet::multinom(target ~ 0 + ., data =  tmp_0_all, trace = FALSE)
  reg_1_sub <- nnet::multinom(target ~ 0 + ., data =  tmp_0_sub, trace = FALSE)

  X_model_matrix_1_all <- stats::model.matrix(reg_1_all)
  xnames_1 <- paste("X", 1:ncol(X_model_matrix_1_all), sep = "")
  znames <- paste("Z", 1:ncol(Z_imp_stand), sep = "")

  #remove unneeded variables
  unneeded <- apply(stats::coefficients(reg_1_sub), 2, function(x) any(is.na(x)))
  xnames_2 <- xnames_1[!unneeded]

  tmp_2_all <- data.frame(target = y_imp)
  tmp_2_all[, xnames_2] <- X_model_matrix_1_all[, !unneeded, drop = FALSE]

  # note: the [] part is only relevant if ncol(X_model_matrix) == 1
  # (see http://stackoverflow.com/questions/40872034/what-happens-when-a-data-frame-gets-new-columns)
  tmp_2_all[, znames] <- Z_imp_stand[, 1:ncol(Z_imp_stand)]
  tmp_2_all[, "ClID"] <- clID

  tmp_2_sub <- tmp_2_all[!missind, , drop = FALSE]


  # -------------- calling the gibbs sampler to get imputation parameters----
  fixformula <- stats::formula(paste("target~ - 1 + ",
                                paste(xnames_2, ":trait", sep = "", collapse = "+"), sep = ""))

  randformula <- stats::formula(paste("~us( - 1 + ", paste(znames, ":trait", sep = "", collapse = "+"),
                                 "):ClID", sep = ""))

  J <- length(table(y_imp)) #number of categories
  number_fix_parameters <- ncol(X_model_matrix_1_all) * (J-1)

  # Get the number of random effects variables
  number_random_effects <- length(znames)
  number_random_parameters <- number_random_effects * (J - 1)
  #Fix residual variance R at 1
  # cf. http://stats.stackexchange.com/questions/32994/what-are-r-structure-g-structure-in-a-glmm


  J_matrix <- array(1, dim = c(J, J) - 1) # matrix of ones
  I_matrix <- diag(J - 1) #identiy matrix

  IJ <- (I_matrix + J_matrix)/J # see Hadfields Course notes p 97
  prior <- list(R = list(V = IJ, fix = TRUE),
                G = list(G1 = list(V = diag(number_random_parameters), nu = 2)))

  MCMCglmm_draws <- MCMCglmm::MCMCglmm(fixed = fixformula,
                                       random = randformula,
                                       rcov = ~us(trait):units,
                                       data = tmp_2_sub,
                                       family = "categorical",
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


  number_of_draws <- nrow(pointdraws)
  select.record <- sample(1:number_of_draws, 1, replace = TRUE)

  # -------------------- drawing samples with the parameters from the gibbs sampler --------
  #now generate new P(Y = A|x * beta) = x*beta/(1+ sum(exp(x*beta))) etc.

  #set up random intercepts and slopes

  y_ret <- data.frame(matrix(nrow = n, ncol = 1))
  ###start imputation

  exp_beta <- array(dim = c(n, J - 1))
  #For each individual (in the rows) a coefficient for the J - 1 categories of the target variable
  #will be saved (in the columns).

  fix_eff_sample <- matrix(xdraws[select.record, ], byrow = TRUE, ncol = J - 1)
  rownames(fix_eff_sample) <- paste("covariate", 1:ncol(X_model_matrix_1_all))
  colnames(fix_eff_sample) <- paste("effect for category", 1:(J-1))
  #xdraws has the following form:
  #c(x1 effect for category 1, x1 effect for cat 2, ..., x1 effect for last cat,
  # x2 effect for cat 1, ..., x2 effect for last cat,
  #... last covariates effect for cat 1, ..., last covariates effect for last cat)

  rand_eff_sample <- matrix(zdraws[select.record, ], nrow = number_clusters)
  rownames(rand_eff_sample) <- paste("cluster", 1:number_clusters)
  
  colnames(rand_eff_sample) <- paste(paste("effect of Z", 1:number_random_effects, 
                                " on category ", sep = ""), 
                                rep(1:(J-1), each = number_random_effects), sep = "")
  
  # zdraws has the following form:
  # effect of Z1 on category 1 in cluster 1, 
  # effect of Z1 on category 1 in cluster 2, 
  # effect of Z1 on category 1 in cluster 3,
  # ... 
  # effect of Z1 on category 2 in cluster 1, 
  # effect of Z1 on category 2 in cluster 2, 
  # effect of Z1 on category 2 in cluster 3, 
  # ...
  # effect of Z2 on category 1 in cluster 1, 
  # effect of Z2 on category 1 in cluster 2, 
  # effect of Z2 on category 1 in cluster 3,
  # ... 
  # effect of Z2 on category 2 in cluster 1, 
  # effect of Z2 on category 2 in cluster 2, 
  # effect of Z2 on category 2 in cluster 3, 
  # ...
  # effect of last covariate on last category (without the reference category) in last cluster.
  
  #rand_eff_sample has then the form
  #row_i is the:
  #effect of Z1 on category 1 in cluster_i, 
  #effect of Z1 on category 2 in cluster_i, 
  #effect of Z2 on category 1 in cluster_i, 
  #effect of Z2 on category 2 in cluster_i.
  
   for(k in 1:ncol(exp_beta)){

    # make for each cluster a matrix with the random effect coefficients
    rand_betas <- rand_eff_sample[, k, drop = FALSE]

    exp_beta[, k] <-
        as.matrix(exp(as.matrix(tmp_2_all[, xnames_2, drop = FALSE]) %*%
                        fix_eff_sample[, k, drop = FALSE] +#fix effects
                        rowSums(tmp_2_all[, znames, drop = FALSE] * rand_betas[clID, , drop = FALSE])))#random effects

      #explanation for the fixed effects part:
      #MCMCglmm_draws$Sol is ordered in the following way: beta_1 for category 1, beta_1 for category_2
      #... beta 1 for category_k, beta_2 for category_1, beta_2 for category_2 ...
      #so we have to skip some values as we proceed by category and not beta.
  }

  y_temp <- array(dim = n)
  for(i in 1:n){

    # ensure that the reference category is in the correct position
    mytmp_prob <- exp_beta[i, ]/(1 + sum(exp_beta[i, ]))
    my_prob <- c(1 - sum(mytmp_prob), mytmp_prob) #the first category is the reference category in MCMCglmm

    y_temp[i] <- sample(levels(y_imp), size = 1, replace = TRUE, prob = my_prob)

  }

  y_ret <- data.frame(y_imp = as.factor(ifelse(is.na(y_imp), y_temp, as.character(y_imp))))



  # --------- returning the imputed data --------------
  return(y_ret)
}


# Generate documentation with devtools::document()
# Build package with devtools::build() and devtools::build(binary = TRUE) for zips
