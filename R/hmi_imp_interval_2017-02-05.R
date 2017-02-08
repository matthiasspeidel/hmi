#' The function to impute interval data variables
#'
#' This functions imputes interval data variables. Those are variables,
#' that consists of a lower and upper (numeric) boundary. Technically
#' those boundaries are cointained in a string, seperated by a semi colon.
#' E.g. if a person reports there income to be something between 3000 and 4000 dollars,
#' its value in the interval covariate would be \code{"3000;4000"}.
#' Left (resp. right) censored data can be donoted by \code{"-Inf;x"} (resp. \code{"x;Inf"}),
#' with \code{x} being the (numeric) observed value.
#' @param y_imp_multi A Vector with the variable to impute.
#' @param X_imp_multi A data.frame with the fixed effects variables.
#' @param M An integer defining the number of imputations that should be made.
#' @return A n x M matrix. Each column is one of M imputed y-variables.
imp_interval <- function(y_imp_multi, X_imp_multi, M){

  ####################################################################
  # BEGIN get starting imputation values by maximizing the likelihood#

  missind <- is.na(y_imp_multi)


  types <- array(dim = ncol(X_imp_multi))
  for(j in 1:length(types)) types[j] <- get_type(X_imp_multi[, j])
  need_stand <- types == "cont"
  categorical <- types == "categorical"

  #remove categories with more than 10 observations as the model in the current form
  #will cause later numerical probles
  too_many_levels <- colnames(X_imp_multi[, categorical, drop = FALSE])[
    apply(X_imp_multi[, categorical, drop = FALSE], 2, function(x) nlevels(factor(x))) > 10]

  X_imp_multi <- X_imp_multi[, !names(X_imp_multi) %in% too_many_levels, drop = FALSE]

  X_imp_multi_stand <- X_imp_multi
  X_imp_multi_stand[, need_stand] <- scale(X_imp_multi[, need_stand])

  # blob has to be numeric, so it must only consists of precise observations
  decomposed <- decompose.interval(interval = y_imp_multi)

  if(any(decomposed$lower > decomposed$upper, na.rm = TRUE)){
    stop("in your interval covariate, some values in the lower bound exceed the upper bound.")
  }

  y_precise_template <- decomposed$precise
  n <- length(y_precise_template)
  #if there are imprecise values only...
  if(all(is.na(y_precise_template))){
    #... the template will be set up with a draw from between the borders
    low_sample <- sample_imp(decomposed$lower)
    up_sample <- sample_imp(decomposed$upper)
    y_precise_template <- msm::rtnorm(n = n, lower = low_sample,
                                      upper = up_sample,
                                      mean = 0,
                                      sd = 1)

      rowMeans(decomposed[, 2:3], na.rm = TRUE)
  }
  blob <- sample_imp(y_precise_template)
  tmp_1 <- data.frame(y = blob)

  # run a linear model to get the suitable model.matrix for imputation of the NAs
  lmstart <- stats::lm(blob ~ 0 + . , data = X_imp_multi_stand)
  X_model_matrix_1 <- stats::model.matrix(lmstart)
  xnames_1 <- paste("X", 1:ncol(X_model_matrix_1), sep = "")
  tmp_1[, xnames_1] <- X_model_matrix_1


  fixformula_1 <- stats::formula(paste("y ~ 0 +", paste(xnames_1, collapse = "+"), sep = ""))


  reg_1 <- stats::lm(fixformula_1, data = tmp_1)

  # remove variables with an NA coefficient

  tmp_2 <- data.frame(y = blob)

  xnames_2 <- xnames_1[!is.na(stats::coefficients(reg_1))]
  tmp_2[, xnames_2] <- X_model_matrix_1[, !is.na(stats::coefficients(reg_1)), drop = FALSE]

  fixformula_2 <- stats::formula(paste("y ~ 0 +", paste(xnames_2, collapse = "+"), sep = ""))
  reg_2 <- stats::lm(fixformula_2, data = tmp_2)
  X_model_matrix_2 <- stats::model.matrix(reg_2)

  max.se <- abs(stats::coef(reg_2) * 3)
  coef.std <- sqrt(diag(stats::vcov(reg_2)))


  includes_unimportants <- any(coef.std > max.se) | any(stats::coef(reg_2) < 1e-03)
  counter <- 0
  while(includes_unimportants & counter <= ncol(X_model_matrix_2)){
    counter <- counter + 1

    X_model_matrix_2 <- as.data.frame(X_model_matrix_2[,
                                      coef.std <= max.se & stats::coef(reg_2) >= 1e-03, drop = FALSE])
    reg_2 <- stats::lm(blob ~ 0 + . , data = X_model_matrix_2)
    #remove regression parameters which have a very high standard error
    max.se <- abs(stats::coef(reg_2) * 3)
    coef.std <- sqrt(diag(stats::vcov(reg_2)))

    includes_unimportants <- any(coef.std > max.se) | any(stats::coef(reg_2) < 1e-03)
  }

  MM_1 <- as.data.frame(X_model_matrix_2)

  # --preparing the ml estimation
  # -define rounding intervals



  #####maximum likelihood estimation using starting values
  ####estimation of the parameters

  #blob <- sample_imp(decompose.interval(y_imp_multi_std)$precise)
  lmstart2 <- stats::lm(blob ~ 0 + ., data = MM_1) # it might be more practical to run the model
  #only based on the observed data, but this could cause some covariates in betastart2 to be dropped
  betastart2 <- as.vector(lmstart2$coef)
  sigmastart2 <- stats::sigma(lmstart2)

  #####maximum likelihood estimation using the starting values

  function_generator <- function(para, X, lower, upper){
    ret <- function(para){
      ret_tmp <- negloglik2_intervalsonly(para = para, X = X,
                            lower = lower, upper = upper)
      return(ret_tmp)
    }
    return(ret)
  }


  #!!! THE STARTING VALUES CAN BE QUITE LOW
  starting_values <- c(betastart2, sigmastart2)

  ###exclude obs below (above) the 0.5% (99.5%) income quantile before maximizing
  ###the likelihood. Reason: Some extrem outliers cause problems during the
  ###maximization

  quants <- stats::quantile(y_precise_template, c(0.005, 0.995), na.rm = TRUE)

  # in X and y_in_negloglik only those observations that are no outliers shall be included.
  # Observations with a missing Y are to be included as well even if they could be an outlier.
  # Therefore w
  keep <- (y_precise_template >= quants[1] & y_precise_template <= quants[2]) |
    is.na(y_precise_template)

  negloglik2_generated <- function_generator(para = starting_values,
                                             X = MM_1[keep, , drop = FALSE],
                                             lower = decomposed$lower[keep],
                                             upper = decomposed$upper[keep])

  m2 <- stats::optim(par = starting_values, negloglik2_generated, method = "BFGS",
                     control = list(maxit = 10000), hessian = TRUE)

  #stats::optim(par = starting_values, negloglik2_generated, method = "Nelder-Mead",
  #             control = list(maxit = 10000), hessian = FALSE)


  par_ml2 <- m2$par
  hess <- m2$hessian

  # link about nearest covariance matrix:
  # http://quant.stackexchange.com/questions/2074/what-is-the-best-way-to-fix-a-covariance-matrix-that-is-not-positive-semi-defi
  # nearPD(hess)$mat
  # isSymmetric(Sigma_ml2)

  Sigma_ml2 <- diag(diag(solve(Matrix::nearPD(hess)$mat)))
  diag(Sigma_ml2) <- pmax(diag(Sigma_ml2), 1e-5)
  # make sure, that the main diagonal elements are non-zero

  ###set starting values equal to the observed income
  ###rounded income will be replaced by imputations later
  #y_std_tmp <- decompose.interval(y_imp_multi_std)$precise

  y_imp <- array(NA, dim = c(n, M))
  imp_tmp <- y_precise_template

  for(j in 1:M){
    ####draw new parameters (because it is a Bayesian imputation)


    pars <- mvtnorm::rmvnorm(1, mean = par_ml2, sigma = Sigma_ml2)
      #first eq on page 63 in Drechsler, Kiesl, Speidel (2015)
      #Can we work with a diagonal matrix as well, or is this too far from the posterior?



    # derive imputation model parameters from previously drawn parameters
    beta_hat <- as.matrix(pars[1:(length(pars) - 1)], ncol = 1)
    sigma_hat <- pars[length(pars)]

    mymean <- as.matrix(MM_1) %*% beta_hat


    #The covariance matrix from equation (3)
    Sigma <- sigma_hat^2


    ###################################
    #BEGIN IMPUTING INTERVALL-DATA AND COMPLETELY MISSING DATA
    #for this purpose we have to replace the lower and upper bounds
    # of those observations with an NA in y_imp_multi by -Inf and Inf

    expanded_lower <- decomposed$lower
    expanded_lower[!is.na(decomposed$precise)] <-
      decomposed$precise[!is.na(decomposed$precise)]
    expanded_lower[is.na(expanded_lower)] <- -Inf

    expanded_upper <- decomposed$upper
    expanded_upper[!is.na(decomposed$precise)] <-
      decomposed$precise[!is.na(decomposed$precise)]
    expanded_upper[is.na(expanded_upper)] <- Inf


    #draw values from the truncated normal distributions
    # the bounds are straight forward for the interval data.
    # for the missing data, the bounds are -Inf and +Inf,
    # which is equivalent to draw from a unbounded normal distribution.
    # for precise observations, the bounds are here set to be NA,
    # resulting in NA draws for those observations.
    # The imputation for precise but rounded data follows in the next section.
    # precise and not rounded data need no impuation at all.
    tnorm_draws <- msm::rtnorm(n = n, lower = expanded_lower,
                                upper = expanded_upper,
                                   mean = as.matrix(MM_1) %*% beta_hat,
                                   sd = sqrt(Sigma))

    y_imp[,j] <- tnorm_draws
  }

  return(y_imp)
}


