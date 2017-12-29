#' The function to impute rounded continuous variables
#'
#' For example the income in surveys is often reported rounded by the respondents.
#' See Drechsler, Kiesl and Speidel (2015) for more details.
#' @param y A Vector with the variable to impute.
#' @param X A data.frame with the fixed effects variables.
#' @param rounding_degrees A numeric vector with the presumed rounding degrees.
#' @references Joerg Drechsler, Hans Kiesl, Matthias Speidel (2015):
#' "MI Double Feature: Multiple Imputation to Address Nonresponse and Rounding Errors in Income Questions".
#' Austrian Journal of Statistics Vol. 44, No. 2, http://dx.doi.org/10.17713/ajs.v44i2.77
#' @return A n x 1 data.frame with the original and imputed values.
imp_roundedcont <- function(y, X, rounding_degrees = c(1, 10, 100, 1000)){

  # ----------------------------- preparing the Y data ------------------
  if(is.factor(y)){
    y <- as.interval(y)
  }

  n <- length(y)

  # ----------------------------- preparing the X data ------------------
  # remove excessive variables
  X1 <- cleanup(X)

  # standardize X
  X1_stand <- stand(X1, rounding_degrees = rounding_degrees)

  #The imputation model of missing values is Y ~ X.
  #In order to get a full model matrix, we need two things
  #1. A place holder ph with an precice structure
  #(meaning that ph is not of class interval. Nevertheless the elements in ph
  #can be an aggregate of imprecise observations (e.g. the mean of lower and upper bound))
  #2. The place holder ph must not contain any NAs, NaNs or Infs.
  decomposed_y <- decompose_interval(interval = y)

  #short check for consistency:
  if(any(decomposed_y[, "lower_general"] > decomposed_y[, "upper_general"], na.rm = TRUE)){
    stop("in your interval covariate, some values in the lower bound exceed the upper bound.")
  }

  # classify the data into the three types of observations:
  # 1. precise data (like 3010 or 3017 - in interval notation "3010;3010", "3017;3017")
  # 2. imprecise data (like "3000;3600")
  # 3. missing data (NA - in interval notation "-Inf;Inf")
  #get the indicator of the missing values
  indicator_precise <- !is.na(decomposed_y[, "precise"])
  indicator_imprecise <- !is.na(decomposed_y[, "lower_imprecise"])
  indicator_missing <- is.infinite(decomposed_y[, "lower_general"]) &
    is.infinite(decomposed_y[, "upper_general"])

  #Preparation for standadizing all observations in y, based on the precise values of y
  #y_precise <- decomposed_y[, "precise"]

  mean_of_y_precise <- mean(decomposed_y[, "precise"], na.rm = TRUE)
  sd_of_y_precise <- stats::sd(decomposed_y[, "precise"], na.rm = TRUE)

  # We intentionally add + 1 because otherwise with the standardized x,
  # the intercept in the regression y ~ x can be exactly 0

  # standardise all observations
  y_stand <- (y - mean_of_y_precise)/sd_of_y_precise + 1

  # standardise the decomposed y
  decomposed_y_stand <- (decomposed_y - mean_of_y_precise)/sd_of_y_precise + 1

  # run a linear model to get the suitable model.matrix for imputation of the NAs
  # Later, another model is run. In many cases, both models are redundant.
  # But in cases with categorical covariates, X_model_matrix_1 will generate
  # additional covariates compared to X_imp_stand.
  # The names of these variables are then stored in tmp_1.
  # Then in the second model it is checked for unneeded variables
  # (e.g. unneeded categories).
  ph_for_y <- sample_imp(rowMeans(decomposed_y_stand[, 4:5]))[, 1]
  df_for_y_on_x <- data.frame(ph_for_y = ph_for_y)

  # run a linear model to get the suitable model.matrix for imputation of the NAs
  xnames_0 <- paste("X", 1:ncol(X1_stand), sep = "")

  df_for_y_on_x[xnames_0] <- X1_stand
  model_y_on_x <-  stats::lm(ph_for_y ~ 0 + . , data = df_for_y_on_x)

  #model matrix
  MM_y_on_x_0 <- stats::model.matrix(model_y_on_x)
  xnames_1 <- paste("X", 1:ncol(MM_y_on_x_0), sep = "")

  df_for_y_on_x <- data.frame(ph_for_y = ph_for_y)
  df_for_y_on_x[, xnames_1] <- MM_y_on_x_0

  safetycounter <- 0
  unneeded <- TRUE
  while(any(unneeded) & safetycounter <= ncol(MM_y_on_x_0)){

    safetycounter <- safetycounter + 1

    # Run another model and...
    reg_1 <- stats::lm(ph_for_y ~ 0 +., data = df_for_y_on_x)
    MM_y_on_x_1 <- stats::model.matrix(reg_1)

    #... remove unneeded variables with an NA coefficient
    unneeded <- is.na(stats::coefficients(reg_1))
    xnames_1 <- colnames(MM_y_on_x_1)[!unneeded]

    df_for_y_on_x <- data.frame(ph_for_y = ph_for_y)
    df_for_y_on_x[, xnames_1] <- MM_y_on_x_1[, !unneeded, drop = FALSE]
  }

  reg_2 <- stats::lm(ph_for_y ~ 0 + ., data = df_for_y_on_x)
  MM_y_on_x_2 <- stats::model.matrix(reg_2)

  # Now check for variables with too much variance
  max.se <- abs(stats::coef(reg_2) * 3)
  coef.std <- sqrt(diag(stats::vcov(reg_2)))

  includes_unimportants <- any(coef.std > max.se)

  safetycounter <- 0
  while(includes_unimportants & safetycounter <= ncol(MM_y_on_x_0)){
    safetycounter <- safetycounter + 1

    xnames_1 <- colnames(MM_y_on_x_2)[coef.std <= max.se]

    df_for_y_on_x <- data.frame(ph_for_y = ph_for_y)
    df_for_y_on_x[, xnames_1] <- MM_y_on_x_2[, xnames_1, drop = FALSE]

    reg_2 <- stats::lm(ph_for_y ~ 0 +., data = df_for_y_on_x)

    MM_y_on_x_2 <- stats::model.matrix(reg_3)
    #check the regression parameters on very high standard errors
    max.se <- abs(stats::coef(reg_2) * 3)
    coef.std <- sqrt(diag(stats::vcov(reg_2)))

    includes_unimportants <- any(coef.std > max.se)
  }


  # --preparing the ml estimation
  # -define rounding intervals

  half_interval_length <- rounding_degrees/2

  # Determine the rounding degrees of the precise observations

  rounding_categories_indicator <- array(dim = c(sum(indicator_precise),
                                                 length(rounding_degrees)))
  for(i in 1:ncol(rounding_categories_indicator)){
    rounding_categories_indicator[, i] <- decomposed_y[indicator_precise, "precise"] %% rounding_degrees[i] == 0
  }

  p <- factor(rowSums(rounding_categories_indicator))

  # Define a matrix for the model p ~ Y + X
  df_for_p_on_y_and_x <- data.frame(ph_for_p = p, df_for_y_on_x[indicator_precise, , drop = FALSE])

  #####maximum likelihood estimation using starting values
  ####estimation of the parameters

  # estimation of the starting values for eta and the thresholds on the x-axis:
  # ordered probit maximum possible rounding on the rounded in income data

  tryCatch(
    {

      #polr throws an warning, if no intercept is included in the model formula
      #(See ?polr)
      #so we add one in the formula and exclude the constant variable in the data.frame
      #before hand.
      constant_variables <- apply(df_for_p_on_y_and_x, 2, function(x) length(unique(x)) == 1)
      df_for_p_on_y_and_x_2 <- df_for_p_on_y_and_x[, !constant_variables, drop = FALSE]
      if(ncol(df_for_p_on_y_and_x_2) == 1){
        probitstart <- MASS::polr("target ~ 0 + .",
                                  data = df_for_p_on_y_and_x,
                                  contrasts = NULL, Hess = TRUE, model = TRUE,
                                  method = "logistic")
      }else{
        probitstart <- MASS::polr("ph_for_p ~ 1 + .",
                                  data = df_for_p_on_y_and_x_2,
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

  gamma1start <- probitstart$coefficients[names(probitstart$coefficients) == "ph_for_y"]

  kstart <- as.vector(probitstart$zeta) # the tresholds (in the summary labeled "Intercepts")
  #explaining the tresholds for the example of rounding degrees 1, 10, 100 and 1000:
  #0 (rounding degree 1), 0|1 (reounding degree 10),  1|2 (100),  2|3 (1000)

  # it might be more practical to run the model
  #only based on the observed data, but this could cause some covariates in betastart2 to be dropped
  betastart <- as.vector(model_y_on_x$coef)
  sigmastart <- sigma(model_y_on_x)

  #####maximum likelihood estimation using the starting values
  #The intercept of the model for y has not be maximized as due to the standardizations
  #of y and x, it's value is exactly 1.
  starting_values <- c(kstart, betastart, gamma1start, sigmastart)

  names(starting_values)[1:length(kstart)] <- paste("threshold", 1:length(kstart), sep = "")
  names(starting_values)[length(kstart) + 1:length(betastart)] <-
    paste("coef_y_on_x", 1:length(betastart), sep = "")
  names(starting_values)[length(kstart) + length(betastart) + 1:length(gamma1start)] <-
    paste("coef_p_on_y_and_x", 1:length(gamma1start), sep = "")
  names(starting_values)[length(starting_values)] <- "sigma"
  ###exclude obs below (above) the 0.5% (99.5%) income quantile before maximizing
  ###the likelihood. Reason: Some extrem outliers cause problems during the
  ###maximization

  quants <- stats::quantile(decomposed_y_stand[indicator_precise, "precise"],
                            c(0.005, 0.995), na.rm = TRUE)

  indicator_outliers <- (decomposed_y_stand[indicator_precise, "precise"] < quants[1] |
                         decomposed_y_stand[indicator_precise, "precise"] > quants[2])

  m2 <- stats::optim(par = starting_values[-(length(kstart) + 1)], negloglik,
                     X_in_negloglik = MM_y_on_x_0,
                     y_precise_stand = decomposed_y_stand[indicator_precise, "precise"],
                     lower_bounds = decomposed_y_stand[indicator_imprecise, 2],
                     upper_bounds = decomposed_y_stand[indicator_imprecise, 3],
                     my_p = as.numeric(as.character(p)),
                     sd_of_y_precise = sd_of_y_precise,
                     rounding_degrees = rounding_degrees,
                     indicator_precise = indicator_precise,
                     indicator_imprecise = indicator_imprecise,
                     indicator_outliers = indicator_outliers,
                     method = "Nelder-Mead",#"BFGS",
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

      Sigma_ml2 <- diag(ncol(hess))
      diag(Sigma_ml2) <- abs(pars)/100

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

 	####draw new parameters (because it is a Bayesian imputation)
  # Boolean value indicating whether the parameters are valid or not
  invalid <- TRUE

  #numerical problems can result in a not positive definite Matrix.
  Sigma_ml3 <- as.matrix(Matrix::nearPD(Sigma_ml2)$mat)
  counter <- 0
  while(invalid & counter < 1000){
    counter <- counter + 1
    pars <- mvtnorm::rmvnorm(1, mean = par_ml2, sigma = Sigma_ml3)
    #first eq on page 63 in Drechsler, Kiesl, Speidel (2015)

    ####test if drawn parameters for the thresholds are in increasing order
    ####and if the standard deviation of the residuals is <= 0
    ####if yes, draw again
    # pars takes the starting values c(kstart, betastart2, gamma1start, sigmastart2)
    invalid <- is.unsorted(pars[1:(length(rounding_degrees) - 1)]) | pars[length(pars)] <= 0
  }

  # derive imputation model parameters from previously drawn parameters
  if(ncol(MM_y_on_x_0) == 1){
    beta_hat <- matrix(1, ncol = 1)
  }else{
    beta_hat <- as.matrix(c(1, pars[length(rounding_degrees):(length(pars) - 2)]), ncol = 1)
  }

  gamma1_hat <- pars[length(pars) - 1]
  sigma_hat <- pars[length(pars)]
  mu_g <- gamma1_hat * (as.matrix(MM_y_on_x_0) %*% beta_hat)
  mu_y <- as.matrix(MM_y_on_x_0) %*% beta_hat

  #The covariance matrix from equation (3)
  Sigma <- matrix(c(1 + gamma1_hat^2 * sigma_hat^2,
                      gamma1_hat * sigma_hat^2, gamma1_hat * sigma_hat^2,
                      sigma_hat^2), nrow = 2)

  ###########################################################
  #BEGIN IMPUTING INTERVAL-DATA AND COMPLETELY MISSING DATA#
  # The imputation for precise but rounded data follows in the next section.
  # precise and not rounded data need no impuation at all.

  lower_general_stand <- decomposed_y_stand[, "lower_general"][indicator_imprecise | indicator_missing]
  upper_general_stand <- decomposed_y_stand[, "upper_general"][indicator_imprecise | indicator_missing]

  #draw values from the truncated normal distributions
  # the bounds are straight forward for the interval data.
  # for the missing data, the bounds are -Inf and +Inf,
  # which is equivalent to draw from a unbounded normal distribution

  mytry_interval <- msm::rtnorm(n = sum(indicator_imprecise | indicator_missing),
                       lower = lower_general_stand,
                       upper = upper_general_stand,
                       mean = mu_y[indicator_imprecise | indicator_missing],
                       sd = sigma_hat)

  # proposed values for imputation
  #do the backtransformation from standardised to unstandardised
  imp_tmp <- decomposed_y[, "precise"]

  imp_tmp[indicator_imprecise | indicator_missing] <-
    (mytry_interval - 1) * sd_of_y_precise + mean_of_y_precise


  ###############################################################################
  ########################### BEGIN UNROUNDING-IMPUTATION########################
  ###define bounds for the rounding basis
  bounds_for_g_hat <- c(-Inf, pars[1:(length(rounding_degrees) - 1)], Inf)
  ###define interval bounds for maximum possible rounding intervals
  #Principally this could be done without standardization, but it makes the following functions
  #work more reliably.
  #If standardization happens, it is important to adjust the parameters accordingly.
  y_lower <- (decomposed_y[indicator_precise, "precise"] -
                half_interval_length[as.numeric(as.character(p))] - mean_of_y_precise)/sd_of_y_precise + 1

  y_upper <- (decomposed_y[indicator_precise, "precise"] +
                half_interval_length[as.numeric(as.character(p))] - mean_of_y_precise)/sd_of_y_precise + 1

  g_upper <- bounds_for_g_hat[as.numeric(as.character(p)) + 1]

  #elements <- cbind(mymean, -Inf, y_lower, g_upper, y_upper)#ORIGINAL
  elements <- cbind(-Inf, mu_g[indicator_precise, 1], g_upper,
                    y_lower, mu_y[indicator_precise, 1], y_upper)

  # Note: we set g_lower to -Inf because we state that a value of 1500 is not necessarily
  # a multiple of 500; it could also be rounded to the next multiple of 10 or even 1.
  colnames(elements) <- c("g_lower", "mean_g","g_upper", "y_lower","mean_y",   "y_upper")

  ###indicator which of the precise observations need to be imputed due to rounding
  #(and not because they are missing)
  rounded <- rep(TRUE, length(p))

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


    #replace the invalid draws with valid ones.
    # For the latent rounding tendency, we use the highest possible rounding tendency
    # For y, we use a uniform sample between the highest and lowest possible
    #bounds of y.
    problematic_draws <- is.na(mytry[, 1])
    problematic_elements <- elements[problematic_draws, , drop = FALSE]

    # check if there are problematic means of g. This is the case if the mean is outside
    # the interval for a possible g.
    toosmall_gs <- problematic_elements[, 2] < problematic_elements[, 1]
    toolarge_gs <- problematic_elements[, 2] > problematic_elements[, 3]

    elements[which(problematic_draws)[toosmall_gs], 2] <-
      elements[which(problematic_draws)[toosmall_gs], 1]

    elements[which(problematic_draws)[toolarge_gs], 2] <-
      elements[which(problematic_draws)[toolarge_gs], 3]


    toosmall_ys <- problematic_elements[, 5] < problematic_elements[, 4]
    toolarge_ys <- problematic_elements[, 5] > problematic_elements[, 6]

    elements[which(problematic_draws)[toosmall_ys], 5] <-
      elements[which(problematic_draws)[toosmall_ys], 4]

    elements[which(problematic_draws)[toolarge_ys], 5] <-
      elements[which(problematic_draws)[toolarge_ys], 6]


    ####get imputed rounding indicator
    round_int <- apply(mytry[, 1, drop = FALSE], 1,
                         function(x) sum(x > bounds_for_g_hat))

    ###get imputed income on original scale
    imp_precise_temp <- (mytry[, 2, drop = FALSE] - 1) * sd_of_y_precise + mean_of_y_precise
    #Store these results as imputation values...
    imp_tmp[indicator_precise][rounded] <- imp_precise_temp

    #... but test if estimated rounding degree and proposed y can explain the observed y.
    # E.g. the estimated rounding degree 10 and the proposed y 2063 doesn't match
    #to an observed value 2100. A degree of 100 would match in this case.
    #If degree and y do match set the value for rounded to FALSE.
    # The remaining (non-matching) observations get a new proposal y and rounding degree.
    domatch <- floor(imp_precise_temp[, 1]/rounding_degrees[round_int] + 0.5) * rounding_degrees[round_int] ==
                          decomposed_y[indicator_precise, "precise"][rounded]
    rounded[rounded][domatch] <- FALSE
  }

  y_ret <- data.frame(y_ret = imp_tmp)

  return(y_ret)
}
