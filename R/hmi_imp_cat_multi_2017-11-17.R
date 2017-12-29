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
#' chains, that are only examined by few iterations (say less than 1000),
#' @param rounding_degrees A numeric vector with the presumed rounding degrees.
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
                      rounding_degrees = c(1, 10, 100, 1000)){

  if(min(table(y_imp)) < 2) {
    stop("Too few observations per category in a categorical target variable.")
  }

  if(!is.factor(y_imp)){
    warning("We suggest to make all your categorical variables to be a factor.")
    y_imp <- as.factor(y_imp)
  }

  # ----------------------------- preparing the X and Z data ------------------

  # remove excessive variables
  X_imp <- cleanup(X_imp)

  # standardise the covariates in X (which are numeric and no intercept)
  X <- stand(X_imp, rounding_degrees = rounding_degrees)

  # -- standardise the covariates in Z (which are numeric and no intercept)
  Z_imp <- cleanup(Z_imp)
  Z <- stand(Z_imp, rounding_degrees = rounding_degrees)

  #the missing indactor indicates, which values of y are missing.
  missind <- is.na(y_imp)
  n <- nrow(X)
  number_clusters <- length(table(clID))


  # ----------------- starting model to get categorical variables as dummies
  # and for eliminating linear dependent variables.
  ph <- sample_imp(y_imp)[, 1]
  YX_0_all <- data.frame(target = ph)
  X2 <- stats::model.matrix(~ 0 + ., data = X)
  xnames_0 <- paste("X", 1:ncol(X2), sep = "")
  YX_0_all[xnames_0] <- X2
  YX_0_sub <- YX_0_all[!missind, , drop = FALSE]

  reg_YX_1_all <- stats::lm(as.numeric(target) ~ 1 +. , data =  YX_0_all)
  X3 <- stats::model.matrix(reg_YX_1_all)[, !is.na(stats::coef(reg_YX_1_all)), drop = FALSE]

  YX_1_all <- data.frame(target = ph)
  xnames_1 <- paste("X", 1:ncol(X3), sep = "")
  YX_1_all[xnames_1] <- X3
  YX_1_sub <- YX_1_all[!missind, , drop = FALSE]


  # -------- better model to remove insignificant variables
  reg_YX_1_all <- nnet::multinom(target ~ 0 + ., data =  YX_1_all, trace = FALSE)
  reg_YX_1_sub <- nnet::multinom(target ~ 0 + ., data =  YX_1_sub, trace = FALSE)

  X_0_all <- stats::model.matrix(reg_YX_1_all)
  # ---------- remove insignificant variables
  # -- check for insignificance
  insignificant <- array(dim = length(reg_YX_1_sub$coefnames))
  for(i in 1:length(insignificant)){
    reg_YX_1_sub_test <- nnet::multinom(paste("target ~ 0 + . -",  reg_YX_1_sub$coefnames[i]),
                                        data =  YX_1_sub, trace = FALSE)
    insignificant[i] <- tryCatch(
      {
      insignificant[i] <- lmtest::lrtest(reg_YX_1_sub, reg_YX_1_sub_test)[[5]][2] > 0.1
      },
      error = function(cond){
        insignificant[i] <- FALSE
      },
      warning = {

      },
      finally = {
      }
    )
  }

  while(any(insignificant)){
    X_0_all <- stats::model.matrix(reg_YX_1_all)[, !insignificant, drop = FALSE]

    YX_1_all <- data.frame(target = ph)
    xnames_1 <- paste("X", 1:ncol(X_0_all), sep = "")
    YX_1_all[xnames_1] <- X_0_all
    YX_1_sub <- YX_1_all[!missind, , drop = FALSE]

    reg_YX_1_all <- nnet::multinom(target ~ 0 + ., data =  YX_1_all, trace = FALSE)
    reg_YX_1_sub <- nnet::multinom(target ~ 0 + ., data =  YX_1_sub, trace = FALSE)

    # -- update insignificance test
    insignificant <- array(dim = length(reg_YX_1_sub$coefnames))
    for(i in 1:length(insignificant)){
      reg_YX_1_sub_test <- nnet::multinom(paste("target ~ 0 + .", "-", reg_YX_1_sub$coefnames[i]),
                                          data =  YX_1_sub,
                                          trace = FALSE)
      insignificant[i] <- tryCatch(
        {
          insignificant[i] <- lmtest::lrtest(reg_YX_1_sub, reg_YX_1_sub_test)[[5]][2] > 0.1
        },
        error = function(cond){
          insignificant[i] <- FALSE
        },
        warning = {

        },
        finally = {
        }
      )
    }

  }

  #Remove correlated variables

  tmptypes <- array(dim = ncol(X_0_all))
  for(i in 1:length(tmptypes)){
    tmptypes[i] <- get_type(X_0_all[, i], rounding_degrees = rounding_degrees)
  }
  somethingcont <- tmptypes %in% c("cont", "binary", "roundedcont", "semicont")
  correlated <- NULL
  if(sum(somethingcont, na.rm = TRUE) >= 2){
    tmp0 <- stats::cor(X_0_all[, somethingcont])
    tmp0[lower.tri(tmp0, diag = TRUE)] <- NA
    tmp1 <- stats::na.omit(as.data.frame(as.table(tmp0)))
    correlated <- unique(as.character(subset(tmp1,
                        abs(tmp1[, 3]) > 0.25 & abs(tmp1[, 3]) < 1)[, 1]))
  }

  if(!is.null(correlated)){
    X_0_all <- X_0_all[, !colnames(X_0_all) %in% correlated]
  }

  YX_1_all <- data.frame(target = ph)
  xnames_1 <- paste("X", 1:ncol(X_0_all), sep = "")
  YX_1_all[xnames_1] <- X_0_all
  YXZ_1_all <- YX_1_all
  znames <- paste("Z", 1:ncol(Z), sep = "")
  YXZ_1_all[znames] <- Z
  YXZ_1_all[, "clID"] <- clID
  YXZ_1_sub <- YXZ_1_all[!missind, , drop = FALSE]

  # -------------- calling the gibbs sampler to get imputation parameters----
  fixformula <- stats::formula(paste("target~ - 1 + ",
                                paste(xnames_1, ":trait", sep = "", collapse = "+"), sep = ""))

  randformula <- stats::formula(paste("~us( - 1 + ", paste(znames, ":trait", sep = "", collapse = "+"),
                                 "):clID", sep = ""))

  J <- length(table(y_imp)) #number of categories
  number_fix_parameters <- ncol(X_0_all) * (J-1)

  # Get the number of random effects variables
  number_random_effects <- length(znames)
  number_random_parameters <- number_random_effects * (J - 1)
  #Fix residual variance R at 1
  # cf. http://stats.stackexchange.com/questions/32994/what-are-r-structure-g-structure-in-a-glmm

  J_matrix <- array(1, dim = c(J, J) - 1) # matrix of ones
  I_matrix <- diag(J - 1) #identiy matrix

  IJ <- (I_matrix + J_matrix)/J # see Hadfield's Course notes p. 97
  prior <- list(R = list(V = IJ, fix = 1),
                G = list(G1 = list(V = diag(number_random_parameters), nu = 2)))

  #Experience showed that categorical models need higher number of Gibbs-Sampler
  #to converge.
  MCMCglmm_draws <- MCMCglmm::MCMCglmm(fixed = fixformula,
                                       random = randformula,
                                       rcov = ~us(trait):units,
                                       data = YXZ_1_sub,
                                       family = "categorical",
                                       verbose = FALSE, pr = TRUE,
                                       prior = prior,
                                       saveX = TRUE, saveZ = TRUE,
                                       nitt = nitt*2,
                                       thin = thin*2,
                                       burnin = burnin*2)

  pointdraws <- MCMCglmm_draws$Sol
  variancedraws <- MCMCglmm_draws$VCV

  xdraws <- pointdraws[, 1:number_fix_parameters, drop = FALSE]
  zdraws <- pointdraws[, number_fix_parameters +
                         1:(number_random_parameters * number_clusters), drop = FALSE]


  number_of_draws <- nrow(pointdraws)
  select.record <- sample(1:number_of_draws, size = 1)

  # -------------------- drawing samples with the parameters from the gibbs sampler --------
  #now generate new P(Y = A|x * beta) = x*beta/(1+ sum(exp(x*beta))) etc.

  #set up random intercepts and slopes

  y_ret <- data.frame(matrix(nrow = n, ncol = 1))
  ###start imputation

  exp_beta <- array(dim = c(n, J - 1))
  #For each individual (in the rows) a coefficient for the J - 1 categories of the target variable
  #will be saved (in the columns).

  fix_eff_sample <- matrix(xdraws[select.record, ], byrow = TRUE, ncol = J - 1)
  rownames(fix_eff_sample) <- paste("covariate", 1:ncol(X_0_all))
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
        as.matrix(exp(as.matrix(YXZ_1_all[, xnames_1, drop = FALSE]) %*%
                        fix_eff_sample[, k, drop = FALSE] +#fix effects
    rowSums(YXZ_1_all[, znames, drop = FALSE] * rand_betas[clID, , drop = FALSE])))#random effects

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

  # --------- returning the imputed data --------------
  ret <- list(y_ret = y_ret, Sol = xdraws, VCV = variancedraws)
  return(ret)
}
