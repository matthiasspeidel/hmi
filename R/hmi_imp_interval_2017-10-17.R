#' The function to impute interval data variables
#'
#' This functions imputes interval data variables. Those are variables,
#' that consists of a lower and upper (numeric) boundary. Technically
#' those boundaries are contained in a string, separated by a semi colon.
#' E.g. if a person reports there income to be something between 3000 and 4000 dollars,
#' its value in the interval covariate would be \code{"3000;4000"}.
#' Left (resp. right) censored data can be denoted by \code{"-Inf;x"} (resp. \code{"x;Inf"}),
#' with \code{x} being the (numeric) observed value.
#' @param y_imp A Vector from the class \code{interval} with the variable to impute.
#' @param X_imp A data.frame with the fixed effects variables.
#' @param rounding_degrees A numeric vector with the presumed rounding degrees.
#' @return A n x 1 data.frame with the original and imputed values.
#' Note that this function won't return \code{interval} data as its purpose is to
#' "break" the interval answers into precise answers.
imp_interval <- function(y_imp, X_imp,
                         rounding_degrees = c(1, 10, 100, 1000)){


  # ----------------------------- preparing the X data ------------------
  # remove excessive variables
  X_imp <- cleanup(X_imp)

  # standardize X
  X_imp_stand <- stand(X_imp, rounding_degrees = rounding_degrees)

  missind <- is.na(y_imp)

  # has to be numeric, so it must only consists of precise observations
  decomposed <- decompose_interval(interval = y_imp)

  if(any(decomposed[, "lower_general"] > decomposed[, "upper_general"], na.rm = TRUE)){
    stop("in your interval covariate, some values in the lower bound exceed the upper bound.")
  }

  y_precise_template <- sample_imp(center.interval(y_imp, inf2NA = TRUE))[, 1]
  n <- length(y_precise_template)

  y_mean <- mean(y_precise_template)
  y_sd <- stats::sd(y_precise_template)
  #add 1 to avoid an exactly zero intercept when modeling y_stand ~ x_stand
  # If both, X and Y, are standardized, the intercept
  #will be exactly 0 and thus not significantly different from 0.
  #So in order to avoid this variable to be removed later in the code, we add +1.
  decomposed_stand <- (decomposed - y_mean)/y_sd + 1
  y_precise_template_stand <- (y_precise_template - y_mean)/y_sd + 1
  #if there are imprecise values only...
  #if(all(is.na(y_precise_template))){
  #... the template will be set up with a draw from between the borders
  low_sample <- decomposed_stand[, "lower_general"]
  up_sample  <- decomposed_stand[, "upper_general"]

  y_precise_template <- msm::rtnorm(n = n, lower = low_sample,
                                      upper = up_sample,
                                      mean = y_precise_template_stand,
                                      sd = 1)
  #}

  ph_stand <- y_precise_template
  tmp_0 <- data.frame(target = ph_stand)

  # run a linear model to get the suitable model.matrix for imputation of the NAs
  xnames_0 <- paste("X", 1:ncol(X_imp_stand), sep = "")
  tmp_0[xnames_0] <- X_imp_stand
  lmstart <- stats::lm(target ~ 0 + . , data = tmp_0)
  X_model_matrix_1 <- stats::model.matrix(lmstart)

  tmp_1 <- data.frame(target = ph_stand)
  xnames_1 <- paste("X", 1:ncol(X_model_matrix_1), sep = "")
  tmp_1[, xnames_1] <- X_model_matrix_1

  reg_1 <- stats::lm(target ~ 0 + ., data = tmp_1)

  # remove variables with an NA coefficient

  tmp_2 <- data.frame(target = ph_stand)

  xnames_2 <- xnames_1[!is.na(stats::coefficients(reg_1))]
  tmp_2[, xnames_2] <- X_model_matrix_1[, !is.na(stats::coefficients(reg_1)), drop = FALSE]

  reg_2 <- stats::lm(target ~ 0 + ., data = tmp_2)
  X_model_matrix_2 <- stats::model.matrix(reg_2)

  max.se <- abs(stats::coef(reg_2) * 3)
  coef.std <- sqrt(diag(stats::vcov(reg_2)))

  includes_unimportants <- any(coef.std > max.se) | any(abs(stats::coef(reg_2)) < 1e-03)
  counter <- 0
  while(includes_unimportants & counter <= ncol(X_model_matrix_2)){
    counter <- counter + 1

    X_model_matrix_2 <- as.data.frame(X_model_matrix_2[,
                                      coef.std <= max.se & stats::coef(reg_2) >= 1e-03, drop = FALSE])
    if(ncol(X_model_matrix_2) == 0){
      reg_2 <- stats::lm(ph_stand ~ 1)
    }else{
      reg_2 <- stats::lm(ph_stand ~ 0 + . , data = X_model_matrix_2)
    }

    #remove regression parameters which have a very high standard error
    max.se <- abs(stats::coef(reg_2) * 3)
    coef.std <- sqrt(diag(stats::vcov(reg_2)))

    includes_unimportants <- any(coef.std > max.se) | any(stats::coef(reg_2) < 1e-03)
  }

  MM_1 <- as.data.frame(X_model_matrix_2)

  tmp_3 <- data.frame(target = ph_stand)

  xnames_3 <- names(MM_1)
  tmp_3[, xnames_3] <- MM_1
  # --preparing the ml estimation
  # -define rounding intervals

  #####maximum likelihood estimation using starting values
  ####estimation of the parameters

  lmstart2 <- stats::lm(target ~ 0 + ., data = tmp_3) # it might be more practical to run the model
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
                                             lower = decomposed_stand[, "lower_imprecise"][keep],
                                             upper = decomposed_stand[, "upper_imprecise"][keep])

  m2 <- stats::optim(par = starting_values, negloglik2_generated,
                     method = "BFGS",#alternative: "Nelder-Mead"
                     control = list(maxit = 10000), hessian = TRUE)


  par_ml2 <- m2$par
  hess <- m2$hessian

  # link about nearest covariance matrix:
  # http://quant.stackexchange.com/questions/2074/what-is-the-best-way-to-fix-a-covariance-matrix-that-is-not-positive-semi-defi
  # nearPD(hess)$mat
  # isSymmetric(Sigma_ml2)


  Sigma_ml2 <- tryCatch(
    {

      Sigma_ml2 <- solve(hess)

    },
    error = function(cond) {
      cat("Hessian matrix couldn't be inverted (in the imputation function of the rounded continuous variable).
              Still, you should get a result, but which needs special attention.\n")
      cat("Here's the original error message:\n")
      cat(as.character(cond))

      Sigma_ml2 <- diag(diag(solve(Matrix::nearPD(hess)$mat)))
      diag(Sigma_ml2) <- pmax(diag(Sigma_ml2), 1e-5)

    },
    warning = function(cond) {
      cat("There seems to be a problem with the Hessian matrix in the imputation of the rounded continuous variable.\n")
      cat("Here's the original warning message:\n")
      cat(as.character(cond))

      Sigma_ml2 <- solve(hess)

    },
    finally = {
    }
  )

  # make sure, that the main diagonal elements are non-zero

  ###set starting values equal to the observed income
  ###rounded income will be replaced by imputations later

  imp_tmp <- y_precise_template


  ####draw new parameters (because it is a Bayesian imputation)

  pars <- mvtnorm::rmvnorm(1, mean = par_ml2, sigma = Sigma_ml2)
  #first eq on page 63 in Drechsler, Kiesl, Speidel (2015)

  # derive imputation model parameters from previously drawn parameters
  beta_hat <- as.matrix(pars[1:(length(pars) - 1)], ncol = 1)
  sigma_hat <- pars[length(pars)]

  mymean <- as.matrix(MM_1) %*% beta_hat


  #The covariance matrix from equation (3)
  Sigma <- sigma_hat^2


  ###################################
  #BEGIN IMPUTING INTERVALL-DATA AND COMPLETELY MISSING DATA
  #for this purpose we have to replace the lower and upper bounds
  # of those observations with an NA in y_imp by -Inf and Inf

  expanded_lower <- decomposed_stand[, "lower_general"]

  expanded_upper <- decomposed_stand[, "upper_general"]


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
                                   mean = mymean,
                                   sd = sqrt(Sigma))

  #undo the standardization
  y_ret <- (tnorm_draws - 1) * y_sd + y_mean

  return(data.frame(y_ret = y_ret))
}


