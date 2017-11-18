#' The function to impute rounded continuous variables
#'
#' For example the income in surveys is often reported rounded by the respondents.
#' See Drechsler, Kiesl and Speidel (2015) for more details.
#' @param y_imp A Vector with the variable to impute.
#' @param X_imp A data.frame with the fixed effects variables.
#' @references Joerg Drechsler, Hans Kiesl, Matthias Speidel (2015):
#' "MI Double Feature: Multiple Imputation to Address Nonresponse and Rounding Errors in Income Questions".
#' Austrian Journal of Statistics Vol. 44, No. 2, http://dx.doi.org/10.17713/ajs.v44i2.77
#' @return A n x 1 data.frame with the original and imputed values.
imp_roundedcont <- function(y_imp, X_imp){

  # ----------------------------- preparing the Y data ------------------
  if(is.factor(y_imp)){
    y_imp <- as.interval(y_imp)
  }

  n <- length(y_imp)

  # ----------------------------- preparing the X data ------------------
  # remove excessive variables
  X_imp <- cleanup(X_imp)

  # standardize X
  X_imp_stand <- stand(X_imp)

  #The imputation model of missing values is Y ~ X.
  #In order to get a full model matrix, we need two things
  #1. A place holder ph with an precice structure
  #(meaning that ph is not of class interval. Nevertheless the elements in ph
  #can be an aggregate of imprecise observations (e.g. the mean of lower and upper bound))
  #2. The place holder ph must not contain any NAs, NaNs or Infs.
  decomposed <- decompose_interval(interval = y_imp)

  #short check for consistency:
  if(any(decomposed[, "lower_general"] > decomposed[, "upper_general"], na.rm = TRUE)){
    stop("in your interval covariate, some values in the lower bound exceed the upper bound.")
  }

  # classify the data into the three types of observations:
  # 1. precise data (like 3010 or 3017 - in interval notation "3010;3010", "3017;3017")
  # 2. imprecise data (like "3000;3600")
  # 3. missing data (NA - in interval notation "-Inf;Inf")
  #get the indicator of the missing values
  indicator_precise <- !is.na(decomposed[, "precise"])
  indicator_imprecise <- !is.na(decomposed[, "lower_imprecise"])
  indicator_missing <- is.infinite(decomposed[, "lower_general"]) &
    is.infinite(decomposed[, "upper_general"])

  # standardise the data
  y_precise <- decomposed[, "precise"]

  mean_y_precise <- mean(y_precise, na.rm = TRUE)
  sd_y_precise <- stats::sd(y_precise, na.rm = TRUE)

  # We intentionally add + 1 because otherwise with the standardized x,
  # the intercept in the regression y ~ x can be exactly 0
  y_imp_precise_stand <- (y_imp - mean_y_precise)/sd_y_precise + 1
  y_precise_stand <- (y_precise - mean_y_precise)/sd_y_precise + 1
  if(is_interval(y_imp_precise_stand)){
    y_imp_precise_stand <- decompose_interval(y_imp_precise_stand)[, "precise"]
  }
  decomposed_stand <- (decomposed - mean_y_precise)/sd_y_precise + 1
  #generate a place holder that will be used in later models to get model.matrices
  ph_stand <- sample_imp(rowMeans(decomposed_stand[, 4:5]))[, 1]

  # run a linear model to get the suitable model.matrix for imputation of the NAs
  # Later, another model is run. In many cases, both models are redundant.
  # But in cases with categorical covariates, X_model_matrix_1 will generate
  # additional covariates compared to X_imp_stand.
  # The names of these variables are then stored in tmp_1.
  # Then in the second model it is checked for unneeded variables
  # (e.g. unneeded categories).
  tmp_0 <- data.frame(target = ph_stand)

  # run a linear model to get the suitable model.matrix for imputation of the NAs
  xnames_0 <- paste("X", 1:ncol(X_imp_stand), sep = "")
  tmp_0[xnames_0] <- X_imp_stand

  reg_0 <- stats::lm(target ~ 0 + . , data = tmp_0)

  # extract the model matrix from the model
  X_model_matrix_1 <- stats::model.matrix(reg_0)
  xnames_1 <- paste("X", 1:ncol(X_model_matrix_1), sep = "")
  tmp_1 <- data.frame(target = ph_stand)
  tmp_1[, xnames_1] <- X_model_matrix_1

  safetycounter <- 0
  unneeded <- TRUE
  while(any(unneeded) & safetycounter <= ncol(X_imp)){

    # Run another model and...
    reg_1 <- stats::lm(target ~ 0 +., data = tmp_1)
    X_model_matrix_1 <- stats::model.matrix(reg_1)

    #... remove unneeded variables with an NA coefficient
    unneeded <- is.na(stats::coefficients(reg_1))
    xnames_1 <- colnames(X_model_matrix_1)[!unneeded]

    tmp_1 <- data.frame(target = ph_stand)
    tmp_1[, xnames_1] <- X_model_matrix_1[, !unneeded, drop = FALSE]

    safetycounter <- safetycounter + 1
  }

  reg_2 <- stats::lm(target ~ 0 + ., data = tmp_1)
  X_model_matrix_2 <- stats::model.matrix(reg_2)

  # Now check for variables with too much variance
  max.se <- abs(stats::coef(reg_2) * 3)
  coef.std <- sqrt(diag(stats::vcov(reg_2)))

  includes_unimportants <- any(coef.std > max.se)

  counter <- 0
  while(includes_unimportants & counter <= ncol(X_model_matrix_2)){
    counter <- counter + 1

    X_model_matrix_2 <- as.data.frame(X_model_matrix_2[, coef.std <= max.se, drop = FALSE])
    lm_less_variables <- stats::lm(ph_stand ~ 0 + . , data = X_model_matrix_2)
    #remove regression parameters which have a very high standard error
    max.se <- abs(stats::coef(lm_less_variables) * 3)
    coef.std <- sqrt(diag(stats::vcov(lm_less_variables)))

    includes_unimportants <- any(coef.std > max.se)
  }

  MM_1 <- as.data.frame(X_model_matrix_2)

  #Define a matrix for the model p ~ Y + X
  MM_p <- cbind(ph_stand, MM_1)


  # --preparing the ml estimation
  # -define rounding intervals

  round_base <- c(1, 10, 100, 1000)
  intervals <- round_base/2

  #check if which observation are rounded
  #Calculate the rounding degree only for those with not an missing value in inc

  #p1 <- y_precise %% 5    ==  0  # divisable by 5
  #p1[is.na(p1)] <- FALSE

  p2 <- y_precise %% 10   ==  0  # divisable by 10
  p2[is.na(p2)] <- FALSE

  #p3 <- y_precise %% 50   ==  0  # etc
  #p3[is.na(p3)] <- FALSE

  p4 <- y_precise %% 100  ==  0  #
  p4[is.na(p4)] <- FALSE

  #p5 <- y_precise %% 500  ==  0  #
  #p5[is.na(p5)] <- FALSE

  p6 <- y_precise %% 1000 ==  0  #
  p6[is.na(p6)] <- FALSE

  p <- factor( p2 +  p4 +  p6,
              levels = c("0", "1", "2", "3"), ordered = TRUE)
  ###indicator which variables need to be imputed because they are rounded (and not because they are missing)
  rounded <- p != 0


  #####maximum likelihood estimation using starting values
  ####estimation of the parameters

  # estimation of the starting values for eta and the thresholds on the x-axis:
  # ordered probit maximum possible rounding on the rounded in income data

  pnames <- colnames(MM_p)
  MM_p1 <- data.frame(target = p)
  MM_p1[, pnames] <- MM_p[, pnames]
  tryCatch(
    {

      #polr throws an warning, if no intercept is included in the model formula
      #(See ?polr)
      #so we add one in the formula and exclude the constant variable in MM_p1
      #before hand.
      #But if now only the target variable is left, because the only "explanatory" covariate
      #was the intercept, we have to use MM_p1 with the intercept in the data
      #but not in the formula as polr would estimate no coefficient.
      constant_variables <- apply(MM_p1, 2, function(x) length(unique(x)) == 1)
      MM_p2 <- MM_p1[indicator_precise, !constant_variables, drop = FALSE]
      if(ncol(MM_p2) == 1){
        probitstart <- MASS::polr("target ~ 0 + .",
                                  data = MM_p1[indicator_precise, , drop = FALSE],
                                  contrasts = NULL, Hess = TRUE, model = TRUE,
                                  method = "probit")
      }else{
        probitstart <- MASS::polr("target ~ 1 + .",
                                  data = MM_p2,
                                  contrasts = NULL, Hess = TRUE, model = TRUE,
                                  method = "probit")
      }


    },
    error = function(cond) {
      cat("We assume that perfect separation occured in your rounded continuous variable, because of too few observations.\n
Consider specifying the variable to be continuous via list_of_types (see ?hmi).\n")
      cat("Here is the original error message:\n")
      cat(as.character(cond))

      return(NULL)
    },
    warning = function(cond) {
      cat("We assume that perfect separation occured in your rounded continuous variable, because of too few observations.\n
Consider specifying the variable to be continuous via list_of_types (see ?hmi).\n")
      cat("Here is the original warning message:\n")
      cat(as.character(cond))

      return(NULL)
    },
    finally = {

    }
  )

  #???More or just one parameter for the rounding degree model???
  gamma1start <- probitstart$coefficients[names(probitstart$coefficients) == "ph_stand"]
  #as.vector(probitstart$coefficients) # the fix effect(s)
  kstart <- as.vector(probitstart$zeta) # the tresholds (in the summary labeled "Intercepts")
  #explaining the tresholds:
  #0 (rounding degree 1), 0|1 (reounding degree 10),  1|2 (100),  2|3 (1000)

  tmp_3 <- data.frame(target = ph_stand)
  xnames_3 <- colnames(MM_1)
  tmp_3[, xnames_3] <- MM_1
  lmstart2 <- stats::lm(target ~ 0 + ., data = tmp_3)
  # it might be more practical to run the model
  #only based on the observed data, but this could cause some covariates in betastart2 to be dropped
  betastart2 <- as.vector(lmstart2$coef)
  sigmastart2 <- summary(lmstart2)$sigma

  #####maximum likelihood estimation using the starting values

  function_generator <- function(para, X, y_in_negloglik, lower, upper,
                                 my_p, mean_y_precise, sd_y_precise){
    ret <- function(para){
      ret_tmp <- negloglik2(para = para, X = X, y_in_negloglik = y_in_negloglik,
                            lower = lower, upper = upper,
                            my_p = my_p,
                            mean_y_precise = mean_y_precise,
                            sd_y_precise = sd_y_precise)
      return(ret_tmp)
    }
    return(ret)
  }


  starting_values <- c(kstart, betastart2, gamma1start, sigmastart2)

  ###exclude obs below (above) the 0.5% (99.5%) income quantile before maximizing
  ###the likelihood. Reason: Some extrem outliers cause problems during the
  ###maximization

  quants <- stats::quantile(y_precise, c(0.005, 0.995), na.rm = TRUE)

  # in X and y_in_negloglik only those observations that are no outliers shall be included.
  # Observations with a missing Y are to be included as well even if they could be an outlier.
  # Therefore w
  keep <- (y_precise>= quants[1] & y_precise <= quants[2]) |
    is.na(y_precise)


  #the interval data have to be standardised as well:
  lower_imprecise_stand <- decomposed_stand[, "lower_imprecise"]
  upper_imprecise_stand <- decomposed_stand[, "upper_imprecise"]

  negloglik2_generated <- function_generator(para = starting_values,
                                             X = MM_1[keep, , drop = FALSE],
                                             y_in_negloglik = y_precise_stand[keep],
                                             lower = lower_imprecise_stand[keep],
                                             upper = upper_imprecise_stand[keep],
                                             my_p = as.numeric(as.character(p[keep])),
                                             mean_y_precise = mean_y_precise,
                                             sd_y_precise = sd_y_precise)

  m2 <- stats::optim(par = starting_values, negloglik2_generated, method = "BFGS",
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
      cat("Here is the original error message:\n")
      cat(as.character(cond))

      Sigma_ml2 <- diag(ncol(hess))/100

    },
    warning = function(cond) {
      cat("There seems to be a problem with the Hessian matrix in the imputation of the rounded continuous variable\n")
      cat("Here is the original warning message:\n")
      cat(as.character(cond))

      Sigma_ml2 <- solve(hess)

    },
    finally = {
    }
  )


  ###set starting values equal to the observed income
  ###rounded income will be replaced by imputations later
 	imp_tmp <- y_precise

 	####draw new parameters (because it is a Bayesian imputation)
  check <- TRUE
  #  counter <- 0
  while(check){
      pars <- mvtnorm::rmvnorm(1, mean = par_ml2, sigma = Sigma_ml2)
      #first eq on page 63 in Drechsler, Kiesl, Speidel (2015)

      ####test if drawn parameters for the thresholds are in increasing order
      ####and if the standard deviation of the residuals is <= 0
      ####if yes, draw again
      # pars takes the starting values c(kstart, betastart2, gammastart, sigmastart2)
      check <- is.unsorted(pars[1:3]) | pars[length(pars)] <= 0
  }


  # derive imputation model parameters from previously drawn parameters
  beta_hat <- as.matrix(pars[4:(length(pars) - 2)], ncol = 1)
  gamma1_hat <- pars[length(pars) - 1]
  sigma_hat <- pars[length(pars)]
  mu_g <- gamma1_hat * as.matrix(MM_1) %*% beta_hat
  mu_y <- as.matrix(MM_1) %*% beta_hat
  mymean <- cbind(mu_g, mu_y)

  #The covariance matrix from equation (3)
  Sigma <- matrix(c(1 + gamma1_hat^2 * sigma_hat^2,
                      gamma1_hat * sigma_hat^2, gamma1_hat * sigma_hat^2,
                      sigma_hat^2), nrow = 2)


  ###########################################################
  #BEGIN IMPUTING INTERVAL-DATA AND COMPLETELY MISSING DATA#
  # The imputation for precise but rounded data follows in the next section.
  # precise and not rounded data need no impuation at all.

  lower_general_stand <- decomposed_stand[, "lower_general"][indicator_imprecise | indicator_missing]
  upper_general_stand <- decomposed_stand[, "upper_general"][indicator_imprecise | indicator_missing]

  #draw values from the truncated normal distributions
  # the bounds are straight forward for the interval data.
  # for the missing data, the bounds are -Inf and +Inf,
  # which is equivalent to draw from a unbounded normal distribution

  mytry <- msm::rtnorm(n = sum(indicator_imprecise | indicator_missing),
                       lower = lower_general_stand,
                        upper = upper_general_stand,
                        mean = as.matrix(MM_1[indicator_imprecise | indicator_missing, , drop = FALSE]) %*%
                         beta_hat,
                        sd = sigma_hat)

  # proposed values for imputation
  #do the backtransformation from standardised to unstandardised
  imp_tmp_imprecise_NA <- (mytry - 1) * sd_y_precise + mean_y_precise

  imp_tmp[indicator_imprecise | indicator_missing] <- imp_tmp_imprecise_NA


  ###############################################################################
  ########################### BEGIN UNROUNDING-IMPUTATION########################
  ###define bounds for the rounding basis
  bounds_hat <- c(-Inf, pars[1:3], Inf)
  ###define interval bounds for maximum possible rounding intervals
  y_lower <- y_precise_stand - intervals[as.numeric(as.character(p)) + 1]/sd_y_precise

  y_upper <- y_precise_stand + intervals[as.numeric(as.character(p)) + 1]/sd_y_precise

  g_lower <- bounds_hat[as.numeric(as.character(p)) + 1]
  g_upper <- bounds_hat[as.numeric(as.character(p)) + 2]

  elements <- cbind(mymean, -Inf, y_lower, g_upper, y_upper)
  #hint: we don't use g_lower because we state that a value of 1500 is not necessarily
  # a multiple of 500; it could also be rounded to the next 10 or oven 1 unit.
  colnames(elements) <- c("mean_g", "mean_y", "g_lower", "y_lower", "g_upper", "y_upper")

  while(any(rounded)){

    ###draw values for g and y from a truncated multivariate normal
    ###drawn y must be between y_lower and y_upper
    ###drawn g must be between g_lower and g_upper
    mytry <- t(apply(elements[rounded, , drop = FALSE],
                    1, sampler, Sigma))

    #It can happen, that rtmvnorm can't sample values from a truncated normal distribution
    #properly. See the following example returning two NaNs
    #instead a values from [0;1]:
    #tmvtnorm::rtmvnorm(1, mean = c(40, 0.5),
    #                   sigma = diag(2),
    #                   lower = c(0, 0),
    #                   upper = c(1, 1),
    #                   algorithm = "gibbs", burn.in.samples = 1000)

    #So, if for individual i no valid value for g or y could be sampled,
    # it could either be because mu_g[i] lies outisde of the interval
    #[g_lower[i];g_upper[i]] or because mu_y[i] outside of y_lower[i];y_upper[i]].
    # We then check whether it is mu_g or mu_y, that loes outside its interval
    #and then replace the corresponding mean
    #by a uniform sample between the lower and the upper bound.


    #replace the draws
    #with valid ones.
    # For the latent rounding tendency, we use the highest possible rounding tendency
    # For y, we use a uniform sample between the highest and lowest possible
    #bounds of y.
    problematic_draws <- is.na(mytry[, 1])
    problematic_elements <- elements[rounded, , drop = FALSE][problematic_draws, , drop = FALSE]

    # check if there are problematic means of g. This is the case if the mean is outside
    # the interval for a possible g.
    toosmall_gs <- problematic_elements[, 1] < problematic_elements[, 3]
    toolarge_gs <- problematic_elements[, 1] > problematic_elements[, 5]

    elements[which(rounded)[which(problematic_draws)[toosmall_gs]], 1] <-
      elements[which(rounded)[which(problematic_draws)[toosmall_gs]], 3]

    elements[which(rounded)[which(problematic_draws)[toolarge_gs]], 1] <-
      elements[which(rounded)[which(problematic_draws)[toolarge_gs]], 5]


    toosmall_ys <- problematic_elements[, 2] < problematic_elements[, 4]
    toolarge_ys <- problematic_elements[, 2] > problematic_elements[, 6]

    elements[which(rounded)[which(problematic_draws)[toosmall_ys]], 2] <-
      elements[which(rounded)[which(problematic_draws)[toosmall_ys]], 4]

    elements[which(rounded)[which(problematic_draws)[toolarge_ys]], 2] <-
      elements[which(rounded)[which(problematic_draws)[toolarge_ys]], 6]


    ####get imputed rounding indicator
    round_int <- apply(mytry[, 1, drop = FALSE], 1,
                         function(x) sum(x > bounds_hat))

    ###get imputed income on original scale
    imp_precise_temp <- (mytry[, 2, drop = FALSE] - 1) * sd_y_precise + mean_y_precise
    #Store these results as imputation values...
    imp_tmp[rounded] <- imp_precise_temp

    #... but test if estimated rounding degree and proposed y can explain the observed y.
    # E.g. the estimated rounding degree 10 and the proposed y 2063 doesn't match
    #to an observed value 2100. A degree of 100 would match in this case.
    #If degree and y do match set the value for rounded to FALSE.
    # The remaining (non-matching) observations get a new proposal y and rounding degree.
    domatch <- floor(imp_precise_temp[, 1]/round_base[round_int] + 0.5) * round_base[round_int] ==
                          y_precise[rounded]
    rounded[as.numeric(names(domatch))[domatch]] <- FALSE


  }

  y_ret <- data.frame(y_ret = imp_tmp)

  return(y_ret)
}


