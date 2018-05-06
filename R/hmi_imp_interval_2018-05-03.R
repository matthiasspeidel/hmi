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
#' @param pvalue A numeric between 0 and 1 denoting the threshold of p-values a variable in the imputation
#' model should not exceed. If they do, they are excluded from the imputation model.
#' @param k An integer defining the allowed maximum of levels in a factor covariate.
#' @return A n x 1 data.frame with the original and imputed values.
#' Note that this function won't return \code{interval} data as its purpose is to
#' "break" the interval answers into precise answers.
imp_interval <- function(y_imp,
                         X_imp,
                         pvalue = 0.2,
                         k = Inf){


  # ----------------------------- preparing the X data ------------------
  # remove excessive variables
  X_imp <- cleanup(X_imp, k = k)

  # standardize X
  X <- stand(X_imp)

  # has to be numeric, so it must only consists of precise observations
  decomposed <- decompose_interval(interval = y_imp)

  # classify the data into the three types of observations:
  # 1. precise data (like 3010 or 3017 - in interval notation "3010;3010", "3017;3017")
  # 2. imprecise data (like "3000;3600")
  # 3. missing data (NA - in interval notation "-Inf;Inf")
  #get the indicator of the missing values
  indicator_precise <- !is.na(decomposed[, "precise"])
  indicator_imprecise <- !is.na(decomposed[, "lower_imprecise"])
  indicator_missing <- is.infinite(decomposed[, "lower_general"]) &
    is.infinite(decomposed[, "upper_general"])

  if(any(decomposed[, "lower_general"] > decomposed[, "upper_general"], na.rm = TRUE)){
    stop("in your interval covariate, some values in the lower bound exceed the upper bound.")
  }

  y_precise_template <- sample_imp(center.interval(y_imp, inf2NA = TRUE))[, 1]
  n <- length(y_precise_template)

  y_mean <- mean(y_precise_template)
  y_sd <- stats::sd(y_precise_template)
  #it can happen that in the template only identical values are present,
  #leading to 0 variance.
  if(y_sd <= 0) y_sd <- 1
  #add 1 to avoid an exactly zero intercept when modeling y_stand ~ x_stand
  # If both, X and Y, are standardized, the intercept
  #will be exactly 0 and thus not significantly different from 0.
  #So in order to avoid this variable to be removed later in the code, we add +1.
  decomposed_stand <- (decomposed - y_mean)/y_sd + 1

  #update y_precise_temple by draws from a truncated normal distribution
  #which ensures, that the template values are in line with the bounds given in y_imp
  y_precise_template <- msm::rtnorm(n = n,
                                    lower = decomposed_stand[, "lower_general"],
                                    upper = decomposed_stand[, "upper_general"],
                                    mean = (y_precise_template - y_mean)/y_sd + 1,
                                    sd = 1)
  #For rare cases where it was not possible to draw a value from the truncated normal distribution,
  #the sample_imp is used.
  y_precise_template <- sample_imp(y_precise_template)[, 1]
  ph <- y_precise_template


  tmp_0_all <- data.frame(target = ph, X)
  xnames_1 <- colnames(X)

  tmp_formula <- paste("target~ 0 + ", paste(xnames_1, collapse = "+"), sep = "")
  reg_1_all <- stats::lm(stats::formula(tmp_formula), data = tmp_0_all)

  X_model_matrix_1_all <- stats::model.matrix(reg_1_all)
  xnames_1 <- paste("X", 1:ncol(X_model_matrix_1_all), sep = "")
  colnames(X_model_matrix_1_all) <- xnames_1

  tmp_0_all <- data.frame(target = ph)
  tmp_0_all[, xnames_1] <- X_model_matrix_1_all

  #From this initial model matrix X_model_matrix_1_all
  #now step by step irrelavant variables are removed.

  # Principally those models are based on precise observations only.
  # But in some data situation, there might be no presice observations, only intervals;
  # then all precise and interval data has to be used.
  # Precise data can be use directly, from imprecise data, a draw from within their bounds is used
  if(sum(indicator_precise) < 30){
    use_indicator <- indicator_precise | indicator_imprecise
  }else{
    use_indicator <- indicator_precise
  }
  X_model_matrix_1_sub <- X_model_matrix_1_all[use_indicator, , drop = FALSE]

  # The first step of the reduction is to remove variables having a non-measurable effect
  # (e.g. due to colinearity) on y.
  # tmp_1 shall include the covariates (like X_model_matrix) and additionally the target variable
  ph_sub <- ph[use_indicator]
  tmp_1_sub <- data.frame(target = ph_sub)
  xnames_1 <- colnames(X_model_matrix_1_sub)
  tmp_1_sub[, xnames_1] <- X_model_matrix_1_sub

  tmp_formula <- paste("target~ 0 + ", paste(xnames_1, collapse = "+"), sep = "")
  reg_1_sub <- stats::lm(stats::formula(tmp_formula) , data = tmp_1_sub)

  #remove unneeded variables
  X_model_matrix_1_sub <- X_model_matrix_1_sub[, !is.na(stats::coefficients(reg_1_sub)),
                                               drop = FALSE]

  # Remove insignificant variables from the imputation model
  check <- TRUE
  while(check){
    tmp_1_sub <- data.frame(target = ph_sub)
    xnames_1 <- colnames(X_model_matrix_1_sub)
    tmp_1_sub[, xnames_1] <- X_model_matrix_1_sub
    tmp_formula <- paste("target~ 0 + ", paste(xnames_1, collapse = "+"), sep = "")
    reg_1_sub <- stats::lm(stats::formula(tmp_formula), data = tmp_1_sub)

    pvalues <- summary(reg_1_sub)$coefficients[, 4]
    insignificant_variables <- which(pvalues > pvalue)
    most_insignificant <- insignificant_variables[which.max(pvalues[insignificant_variables])]

    if(length(most_insignificant) == 0){
      check <- FALSE
    }else{
      #.. drop the insignificant variable from the model.matrix, but only if at least 1 variable remains
      tmp_MM <- stats::model.matrix(reg_1_sub)[, -most_insignificant, drop = FALSE]
      if(ncol(tmp_MM) == 0){
        check <- FALSE
      }else{
        X_model_matrix_1_sub <- tmp_MM
      }
    }
  }

  tmp_2_all <- tmp_0_all[, colnames(tmp_1_sub), drop = FALSE]

  #####maximum likelihood estimation using starting values
  ####estimation of the parameters

  lmstart2 <- stats::lm(target ~ 0 + ., data = tmp_1_sub)
  betastart2 <- as.vector(lmstart2$coef)
  sigmastart2 <- stats::sigma(lmstart2)

  #####maximum likelihood estimation using the starting values

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

  m2 <- stats::nlm(f = negloglik2_intervalsonly, p = starting_values,
                   parnames = names(starting_values),
                   X = as.matrix(tmp_2_all[, xnames_1, drop = FALSE])[keep, , drop = FALSE],
                   lower_bounds = decomposed_stand[, "lower_imprecise"][keep],
                   upper_bounds = decomposed_stand[, "upper_imprecise"][keep],
                   hessian = TRUE, gradtol = 1e-4, steptol = 1e-4)


  par_ml2 <- m2$estimate
  names(par_ml2) <- names(starting_values)
  hess <- m2$hessian

  if(m2$code > 2){
    warning(paste("Likelihood-optimization for rounded continuous variable", colnames(y_df), "failed."))
  }

  # link about nearest covariance matrix:
  # http://quant.stackexchange.com/questions/2074/what-is-the-best-way-to-fix-a-covariance-matrix-that-is-not-positive-semi-defi
  # nearPD(hess)$mat
  # isSymmetric(Sigma_ml2)


  Sigma_ml2 <- tryCatch(
    {
      solve(hess)
    },
    error = function(cond) {
      cat("Hessian matrix couldn't be inverted (in the imputation function for interval variables).
          Still, you should get a result, but which needs special attention.\n")

      tmp <- matrix(0, nrow = length(par_ml2), ncol = length(par_ml2))
      diag(tmp) <- abs(par_ml2)/100
      return(tmp)

    },
    warning = function(cond) {
      cat("There seems to be a problem with the Hessian matrix in the imputation of the rounded continuous variable\n")
      cat("Here is the original warning message:\n")
      cat(as.character(cond))
      return(solve(hess))
    },
    finally = {
    }
  )

  # make sure, that the main diagonal elements are non-zero

  ###set starting values equal to the observed income
  ###rounded income will be replaced by imputations later

  imp_tmp <- y_precise_template


  ####draw new parameters (because it is a Bayesian imputation)
  #a negative draw for the variance parameter has to be rejected
  Sigma_ml3 <- as.matrix(Matrix::nearPD(Sigma_ml2)$mat)
  invalid <- TRUE
  counter <- 0
  while(invalid & counter < 1000){
    counter <- counter + 1
    pars <- mvtnorm::rmvnorm(1, mean = par_ml2, sigma = Sigma_ml3)
    invalid <- pars[length(pars)] <= 0
  }

  #first eq on page 63 in Drechsler, Kiesl, Speidel (2015)

  # derive imputation model parameters from previously drawn parameters
  #The covariance matrix from equation (3)
  beta_hat <- as.matrix(pars[1:(length(pars) - 1)], ncol = 1)
  Sigma <- pars[length(pars)]^2

  mymean <- as.matrix(tmp_2_all[, xnames_1, drop = FALSE]) %*% beta_hat

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
  tnorm_draws <- msm::rtnorm(n = n,
                             lower = expanded_lower,
                             upper = expanded_upper,
                             mean = mymean,
                             sd = sqrt(Sigma))

  #undo the standardization
  y_ret <- (tnorm_draws - 1) * y_sd + y_mean

  return(data.frame(y_ret = y_ret))
}


