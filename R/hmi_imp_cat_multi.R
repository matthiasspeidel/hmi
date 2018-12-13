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
imp_cat_multi <- function(y_imp,
                      X_imp,
                      Z_imp,
                      clID,
                      nitt = 22000,
                      burnin = 2000,
                      thin = 20,
                      pvalue = 0.2,
                      k = Inf){

  if(min(table(y_imp)) < 2) {
    stop("Too few observations per category in a categorical target variable.")
  }

  orgclass <- class(y_imp)
  if(!is.factor(y_imp)){
    warning("We suggest to make all your categorical variables to be a factor.")
    y_imp <- as.factor(y_imp)
  }

  # ----------------------------- preparing the X and Z data ------------------

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

  YZ <- data.frame(target = stats::rnorm(n), Z)
  #remove intercept variable
  YZ <- YZ[, apply(YZ, 2, get_type) != "intercept", drop = FALSE]

  Z2 <- stats::model.matrix(stats::lm("target ~ 1 + .", data = YZ))

  tmp_0_all <- data.frame(target = ph)
  xnames_1 <- paste("X", 1:ncol(X), sep = "")
  znames_1 <- paste("Z", 1:ncol(Z2), sep = "")
  colnames(Z2) <- znames_1

  tmp_0_all[, xnames_1] <- X
  tmp_0_all[, znames_1] <- Z2
  tmp_0_all[, "clID"] <- clID

  tmp_formula <- paste("target ~ 0 +",
                       paste(xnames_1, collapse = "+"))

  # If both, an intercept variable and a categorical variable are present in the data,
  # One variable in the model is redundant. This is handled later in the code, so here
  # the default message from lmer is bothering and therefore suppressed.
  oldw <- getOption("warn")
  options(warn = -1)
  suppressMessages(reg_1_all <- nnet::multinom(stats::formula(tmp_formula), data = tmp_0_all,
                                               trace = FALSE))
  options(warn = oldw)

  X_model_matrix_1_all <- stats::model.matrix(reg_1_all)
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
                       paste(xnames_1, collapse = "+"))

  oldw <- getOption("warn")
  options(warn = -1)
  suppressMessages(reg_1_sub <- nnet::multinom(stats::formula(tmp_formula), data = tmp_1_sub,
                                            trace = FALSE))
  options(warn = oldw)

  #remove unneeded variables
  tmp <- stats::coefficients(reg_1_sub)
  X_model_matrix_1_sub <- X_model_matrix_1_sub[, !apply(tmp, 2, function(x) any(is.na(x))),
                                               drop = FALSE]

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
                         paste(xnames_1, collapse = "+"))

    oldw <- getOption("warn")
    options(warn = -1)
    suppressMessages(reg_1_sub <- nnet::multinom(stats::formula(tmp_formula), data = tmp_1_sub,
                                                 trace = FALSE))
    options(warn = oldw)
    z <- summary(reg_1_sub)$coefficients / summary(reg_1_sub)$standard.errors
    pvalues <- apply((1 - stats::pnorm(abs(z)))*2, 2, min)
    insignificant_variables <- which(pvalues > pvalue)
    most_insignificant <- insignificant_variables[which.max(pvalues[insignificant_variables])]

    if(length(most_insignificant) == 0){
      check <- FALSE
    }else{
      tmp <- stats::model.matrix(reg_1_sub) #if an additional intercept variable is included by the model
      #we cannot run stats::model.matrix(reg_1_sub)[, -most_insignificant]
      #Because most_insignificant refers to a situation without an intercept variable.

      #.. drop the insignificant variable from the model.matrix, but only if at least 1 variable remains
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

  if(length(xnames_1) == 0){
    fixformula <- stats::formula("target ~ 1:trait")
  }else{
    fixformula <- stats::formula(paste("target~ - 1 + ",
                                       paste(xnames_1, ":trait", sep = "", collapse = "+"), sep = ""))
  }

  if(length(znames_1) == 0){
    randformula <- stats::formula("~us(1:trait):clID")
  }else{
    randformula <- stats::formula(paste("~us( - 1 + ", paste(znames_1, ":trait", sep = "", collapse = "+"),
                                        "):clID", sep = ""))
  }

  # -------------- calling the gibbs sampler to get imputation parameters----
  J <- length(unique(ph_sub)) #number of categories
  number_fix_parameters <- ncol(X_model_matrix_1_sub) * (J-1)

  # Get the number of random effects variables
  number_random_effects <- length(znames_1)
  number_random_parameters <- number_random_effects * (J - 1)
  #Fix residual variance R at 1
  # cf. http://stats.stackexchange.com/questions/32994/what-are-r-structure-g-structure-in-a-glmm

  J_matrix <- array(1, dim = c(J, J) - 1) # matrix of ones
  I_matrix <- diag(J - 1) #identiy matrix

  IJ <- (I_matrix + J_matrix)/J # see Hadfield's Course notes p. 97
  prior <- list(R = list(V = IJ, fix = 1),
                G = list(G1 = list(V = diag(number_random_parameters), nu = J+number_random_effects)))
  # Note: the success of MCMCglmm highly depends on the size of the parameter nu,
  # specifying the degrees of freedom of the inverse Wishart-Distribution


  #Experience showed that categorical models need higher number of Gibbs-Sampler
  #to converge.
  MCMCglmm_draws <- MCMCglmm::MCMCglmm(fixed = fixformula,
                                       random = randformula,
                                       rcov = ~us(trait):units,
                                       data = YXZ_2_sub,
                                       family = "categorical",
                                       verbose = FALSE, pr = TRUE,
                                       prior = prior,
                                       saveX = TRUE, saveZ = TRUE,
                                       nitt = nitt * 2,
                                       thin = thin * 2,
                                       burnin = burnin * 2)

  pointdraws <- MCMCglmm_draws$Sol
  variancedraws <- MCMCglmm_draws$VCV

  number_fix_parameters <-  ncol(X_model_matrix_1_sub)*(J-1)
  xdraws <- pointdraws[, 1:number_fix_parameters, drop = FALSE]
  zdraws <- pointdraws[, number_fix_parameters +
                         1:(number_random_parameters * number_clusters), drop = FALSE]

  #xdraws has the following form (column wise):
  #c(x1 effect for category 1, x1 effect for cat 2, ..., x1 effect for last cat,
  # x2 effect for cat 1, ..., x2 effect for last cat,
  #... last covariates effect for cat 1, ..., last covariates effect for last cat)

  # zdraws has the following form (column wise):
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

  #old: number_random_parameters
  number_of_draws <- nrow(pointdraws)
  select.record <- sample(1:number_of_draws, size = 1)

  # -------------------- drawing samples with the parameters from the gibbs sampler --------
  #now generate new P(Y = A|x * beta) = x*beta/(1+ sum(exp(x*beta))) etc.

  #set up random intercepts and slopes

  y_ret <- data.frame(matrix(nrow = n, ncol = 1))
  ###start imputation

  # Set up a matrix including the chances of landing in category tau, compared to the reference category.
  # Each row is for one individual. The columns is for the chances of category tau compared to the reference category.
  # Again in other words: For each individual (in the rows) a coefficient for the J - 1 categories of the target variable
  # will be saved (in the columns).
  exp_beta <- array(dim = c(n, J - 1))


  fix_eff_sample <- matrix(xdraws[select.record, , drop = FALSE], byrow = TRUE, ncol = J - 1)
  rownames(fix_eff_sample) <- paste("covariate", 1:ncol(X_model_matrix_1_sub))
  colnames(fix_eff_sample) <- paste("effect for category", 1:(J-1))


  rand_eff_sample <- matrix(zdraws[select.record, , drop = FALSE], nrow = number_clusters)
  rownames(rand_eff_sample) <- paste("cluster", 1:number_clusters)
  colnames(rand_eff_sample) <- paste("dummy", 1:ncol(rand_eff_sample))
  counter <- 0
  for(zindex in 1:number_random_effects){
    for(catindex in 1:(J-1)){
      counter <- counter + 1
      colnames(rand_eff_sample)[counter] <- paste("effect of Z", zindex, " on category ", catindex, sep = "")
    }
  }

  #rand_eff_sample has then the form
  #row_i is the:
  #effect of Z1 on category 1 in cluster_i
  #effect of Z1 on category 2 in cluster_i
  #...
  #effect of Z2 on category 1 in cluster_i
  #effect of Z2 on category 2 in cluster_i
  #...
  #effect of last variable on category 1 in cluster_i
  #effect of last variable on category 2 in cluster_i
  #...
  #effect of last variable on last category in cluster_i
  # Note: by "last category" the J-1 th category is ment (we do not consider the reference category)

   for(K in 1:ncol(exp_beta)){

    # make for each cluster a matrix with the random effect coefficients
    rand_betas <- rand_eff_sample[, grep(paste("category", K), colnames(rand_eff_sample)), drop = FALSE]

    exp_beta[, K] <-
        as.matrix(exp(as.matrix(tmp_2_all[, xnames_1, drop = FALSE]) %*%
                        fix_eff_sample[, K, drop = FALSE] + # fix effects
    rowSums(tmp_2_all[, znames_1, drop = FALSE] * rand_betas[clID, , drop = FALSE])))# random effects

      #explanation for the fixed effects part:
      #MCMCglmm_draws$Sol is ordered in the following way: beta_1 for category 1, beta_1 for category_2
      #... beta 1 for category_k, beta_2 for category_1, beta_2 for category_2 ...
      #so we have to skip some values as we proceed by category and not beta.
  }

  y_temp <- array(dim = n)

  for(i in 1:n){

    # ensure that the reference category is in the correct position
    mytmp_prob <- exp_beta[i, ]/(1 + sum(exp_beta[i, ]))
    my_prob <- c(max(0, 1 - sum(mytmp_prob)), mytmp_prob) #the first category is the reference category in MCMCglmm

    y_temp[i] <- sample(levels(y_imp), size = 1, prob = my_prob)

  }

  y_ret <- data.frame(y_ret = as.factor(ifelse(is.na(y_imp), y_temp, as.character(y_imp))))
  if(orgclass == "character") y_ret$y_ret <- as.character(y_ret$y_ret)

  # --------- returning the imputed data --------------
  ret <- list(y_ret = y_ret, Sol = xdraws, VCV = variancedraws)
  return(ret)
}
