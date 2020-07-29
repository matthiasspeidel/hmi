#' The function to impute rounded continuous variables
#'
#' For example the income in surveys is often reported rounded by the respondents.
#' See Drechsler, Kiesl and Speidel (2015) for more details.
#' @param y_df A data.frame with the variable to impute.
#' @param X_imp A data.frame with the fixed effects variables explaining y_df.
#' @param PSI A data.frame with the variables explaining the latent rounding tendency G.
#' @param pvalue A numeric between 0 and 1 denoting the threshold of p-values a variable in the imputation
#' model should not exceed. If they do, they are excluded from the imputation model.
#' @param k An integer defining the allowed maximum of levels in a factor covariate.
#' @param rounding_degrees A numeric vector with the presumed rounding degrees for Y.
#' @references Joerg Drechsler, Hans Kiesl, Matthias Speidel (2015):
#' "MI Double Feature: Multiple Imputation to Address Nonresponse and Rounding Errors in Income Questions".
#' Austrian Journal of Statistics Vol. 44, No. 2, http://dx.doi.org/10.17713/ajs.v44i2.77
#' @return A n x 1 data.frame with the original and imputed values.
imp_roundedcont <- function(y_df,
                            X_imp,
                            PSI,
                            pvalue = 0.2,
                            k = Inf,
                            rounding_degrees = NULL){

  # ----------------------------- preparing the Y data ------------------
  y <- y_df[, 1]
  if(is.factor(y)){
    y <- as.interval(y)
  }

  #currently the method does not consider 0 (no rounding)
  no_rounding_at_all <- y %% 1 != 0
  y <- round(y)
  n <- length(y)

  # ----------------------------- preparing the X data ------------------
  # remove excessive variables
  X_imp <- cleanup(X_imp, k = k)

  # standardize X
  X <- stand(X_imp)

  #If no default rounding_degrees were given, suggest_rounding_degrees suggets them
  if(is.null(rounding_degrees)){
    rounding_degrees <- suggest_rounding_degrees(y)
  }

  #... if they are still NULL, then c(1, 10, 100, 1000) is used as a default
  if(is.null(rounding_degrees)){
    rounding_degrees <- c(1, 10, 100, 1000)
  }

  #The imputation model of missing values is Y ~ X.
  #In order to get a full model matrix, we need two things
  #1. A place holder ph with a precice structure
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

  # Standardise all observations:
  # We intentionally add + 1 because otherwise with the standardized x,
  # the intercept in the regression y ~ x can be exactly 0.
  y_stand <- (y - mean_of_y_precise)/sd_of_y_precise + 1
  y_precise_template <- sample_imp(center_interval(y, inf2NA = TRUE))[, 1]
  y_precise_template_stand <- (y_precise_template - mean_of_y_precise)/sd_of_y_precise + 1
  # standardise the decomposed y
  decomposed_y_stand <- (decomposed_y - mean_of_y_precise)/sd_of_y_precise + 1


  low_sample <- decomposed_y_stand[, "lower_general"]
  up_sample  <- decomposed_y_stand[, "upper_general"]

  y_precise_template <- msm::rtnorm(n = n, lower = low_sample,
                                    upper = up_sample,
                                    mean = y_precise_template_stand,
                                    sd = 1)

  ph <- y_precise_template

  # prepare the estimation of Y ~ X
  Y_X_0_all <- data.frame(target = ph, X)
  xnames_1 <- colnames(X)

  Y_X_formula <- paste("target~ 0 + ", paste(xnames_1, collapse = "+"), sep = "")
  reg_Y_X_1_all <- stats::lm(stats::formula(Y_X_formula), data = Y_X_0_all)

  X_model_matrix_1_all <- stats::model.matrix(reg_Y_X_1_all)
  xnames_1 <- paste("X", 1:ncol(X_model_matrix_1_all), sep = "")
  colnames(X_model_matrix_1_all) <- xnames_1

  Y_X_0_all <- data.frame(target = ph)
  Y_X_0_all[, xnames_1] <- X_model_matrix_1_all

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
  X_Y_1_sub <- data.frame(target = ph_sub)
  xnames_1 <- colnames(X_model_matrix_1_sub)
  X_Y_1_sub[, xnames_1] <- X_model_matrix_1_sub

  Y_X_formula <- paste("target~ 0 + ", paste(xnames_1, collapse = "+"), sep = "")
  reg_Y_X_1_sub <- stats::lm(stats::formula(Y_X_formula) , data = X_Y_1_sub)

  #remove unneeded variables
  X_model_matrix_1_sub <- X_model_matrix_1_sub[, !is.na(stats::coefficients(reg_Y_X_1_sub)),
                                               drop = FALSE]

  # Remove insignificant variables from the imputation model
  check <- TRUE
  while(check){
    X_Y_1_sub <- data.frame(target = ph_sub)
    xnames_1 <- colnames(X_model_matrix_1_sub)
    X_Y_1_sub[, xnames_1] <- X_model_matrix_1_sub
    Y_X_formula <- paste("target~ 0 + ", paste(xnames_1, collapse = "+"), sep = "")
    reg_Y_X_1_sub <- stats::lm(stats::formula(Y_X_formula), data = X_Y_1_sub)

    pvalues <- summary(reg_Y_X_1_sub)$coefficients[, 4]
    insignificant_variables <- which(pvalues > pvalue)
    most_insignificant <- insignificant_variables[which.max(pvalues[insignificant_variables])]

    if(length(most_insignificant) == 0){
      check <- FALSE
    }else{
      #.. drop the insignificant variable from the model.matrix, but only if at least 1 variable remains
      tmp_MM <- stats::model.matrix(reg_Y_X_1_sub)[, -most_insignificant, drop = FALSE]
      if(ncol(tmp_MM) == 0){
        check <- FALSE
      }else{
        X_model_matrix_1_sub <- tmp_MM
      }
    }
  }

  Y_X_2_all <- Y_X_0_all[, colnames(X_Y_1_sub), drop = FALSE]


  betastart <- as.vector(reg_Y_X_1_sub$coef)
  sigmastart <- stats::sigma(reg_Y_X_1_sub)
  if(is.na(sigmastart)) sigmastart <- 1

  ##### Preparation of rounding degrees start model #######

  half_interval_length <- rounding_degrees/2

  # Determine the rounding degrees of the precise observations

  rounding_categories_indicator <- array(0, dim = sum(indicator_precise))


  for(i in 1:length(rounding_degrees)){
    rounding_categories_indicator <- ifelse(decomposed_y[indicator_precise, "precise"] %% rounding_degrees[i] == 0,
                                            i, rounding_categories_indicator)
  }

  g <- factor(rounding_categories_indicator, ordered = TRUE)

  ##### preparing the PSI data (the variables explaining the rounding tendency G) ####
  #check if y is part of PSI
  is_y_in_PSI <- colnames(y_df)[1] %in% colnames(PSI)
  # reason: it is a special covariate in PSI.
  # It has a separat coefficient and is the only variable allowed being an interval variable.

  PSIvars_being_intervals <- apply(PSI, 2, is_interval)
  PSI <- PSI[, !PSIvars_being_intervals, drop = FALSE]

  # This dataset must not include any NA, except y can include missing values, as
  # the rounding tendency is only calculated based on precise observations
  PSIvars_having_NAs <- apply(PSI, 2, function(x) any(is.na(x)))
  PSI <- PSI[, !PSIvars_having_NAs, drop = FALSE]
  #if y was a part of PSI, its precise part is included
  if(is_y_in_PSI){
    PSI[, colnames(y_df)[1]] <- decomposed_y[, "precise", drop = FALSE]
  }
  #standardize PSI
  PSI <- stand(PSI)

  # Define a matrix for the model G ~ PSI
  df_for_g_sub <- data.frame(target = g, PSI[indicator_precise, , drop = FALSE])
  #note: even if Y is not part of PSI, the conditioning on precise y is needed:
  #imprecise Y (like "2000;3000") don't have a rounding degree.

  #####maximum likelihood estimation using starting values
  ####estimation of the parameters

  # estimation of the starting values for eta and the thresholds on the x-axis:
  # ordered probit maximum possible rounding on the rounded in income data

  probitstart <- tryCatch(
    {

      #polr throws an warning, if no intercept is included in the model formula
      #(See ?polr)
      #so we add one in the formula and exclude the constant variable in the data.frame
      #before hand.
      constant_variables <- apply(df_for_g_sub, 2, function(x) length(unique(x)) == 1)
      df_for_g_sub2 <- df_for_g_sub[, !constant_variables, drop = FALSE]
      if(ncol(df_for_g_sub2) == 1){  # only target variable is left
        probitstart <- ordinal::clm("target ~ 1", data = df_for_g_sub2)
      }else{
        probitstart <- ordinal::clm("target ~ 1 + .", data = df_for_g_sub2)

      }
      probitstart

    },
    error = function(cond) {
      stop("We assume that perfect separation occured in your rounded continuous variable, because of too few observations.\n
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
      constant_variables <- apply(df_for_g_sub, 2, function(x) length(unique(x)) == 1)
      df_for_g_sub2 <- df_for_g_sub[, !constant_variables, drop = FALSE]
      if(ncol(df_for_g_sub2) == 1){
        probitstart <- ordinal::clm("target ~ 0 + .", data = df_for_g_sub)
      }else{
        probitstart <-  ordinal::clm("target ~ 1 + .", data = df_for_g_sub2)

      }
      return(probitstart)
    },
    finally = {

    }
  )


  PSI_as_MM_sub <- stats::model.matrix(probitstart)[[1]]

  #remove the intercept variables
  PSI_as_MM_sub <- PSI_as_MM_sub[, -grep("(Intercept)", colnames(PSI_as_MM_sub)), drop = FALSE]

  #starting values for the thresholds:
  tmp <- stats::coef(probitstart)
  kstart <- as.vector(tmp[grep("[|]", names(tmp))])
  # the tresholds (in the summary labeled "Intercepts")
  #explaining the tresholds for the example of rounding degrees 1, 10, 100 and 1000:
  #0 (rounding degree 1), 0|1 (reounding degree 10),  1|2 (100),  2|3 (1000)

  # other regression parameters in the rounding model
  gammastart <- -tmp[-grep("[|]", names(tmp))] # note the minus since clm uses this parametrization


  gamma1start <- NULL
  gamma1name <- NULL

  # If Y is part of PSI, its coefficient gamma1 has to be treated seperately
  if(is_y_in_PSI){
    gamma1start <- gammastart[names(gammastart) == colnames(y_df)[1]]
    gamma1name <- "gamma1"
  }


  # (For example, gamma1 is part of the covariance matrix in equation 3 in Drechsler, Kiesl, Speidel, 2015)
  gammastart_without_y <- gammastart[names(gammastart) != colnames(y_df)[1]]
  gammastart_without_y_name <- paste("coef_g_on_psi", 1:length(gammastart_without_y), sep = "")
  #If Y would be the only variable, set gamma_without_y_name NULL
  if(length(gammastart_without_y) == 0){
    gammastart_without_y <- NULL
    gammastart_without_y_name <- NULL
  }

  #####maximum likelihood estimation using the starting values
  #The intercept of the model for y has not be maximized as due to the standardizations
  #of y and x, it's value is exactly 1. later it will be removed.
  starting_values <- c(kstart, betastart, gamma1start, gammastart_without_y, sigmastart)
  #tau is not included as it has to be fixed at 1 to make the ordered probit model identifiable
  #c.f. p.62 in Drechsler, Kiesl, Speidel (2015)


  names(starting_values)[1:length(kstart)] <- paste("threshold", 1:length(kstart), sep = "")

  names(starting_values)[length(kstart) + 1:length(betastart)] <-
    paste("coef_y_on_x", 1:length(betastart), sep = "")

  if(length(c(gamma1name, gammastart_without_y_name)) > 0){
    names(starting_values)[length(kstart) + length(betastart) +
                             1:length(gammastart)] <- c(gamma1name, gammastart_without_y_name)
  }


  names(starting_values)[length(starting_values)] <- "sigma"
  ###exclude obs below (above) the 0.5% (99.5%) income quantile before maximizing
  ###the likelihood. Reason: Some extrem outliers probably cause problems during the
  ###maximization

  quants <- stats::quantile(decomposed_y_stand[indicator_precise, "precise"],
                            c(0.005, 0.995), na.rm = TRUE)

  indicator_outliers <- (decomposed_y_stand[indicator_precise, "precise"] < quants[1] |
                         decomposed_y_stand[indicator_precise, "precise"] > quants[2])

  PSI_in_negloglik <- PSI_as_MM_sub[, colnames(PSI_as_MM_sub) != colnames(y_df)[1], drop = FALSE]
  oldw <- getOption("warn")
  options(warn = -1)
  on.exit(options(warn = oldw))

  m2 <- stats::nlm(f = negloglik, p = starting_values,
                   parnames = names(starting_values),
                   X_in_negloglik = Y_X_2_all[ , xnames_1, drop = FALSE],
                   PSI_in_negloglik = PSI_in_negloglik,
                   y_precise_stand = decomposed_y_stand[indicator_precise, "precise"],
                   lower_bounds = decomposed_y_stand[indicator_imprecise, 2],
                   upper_bounds = decomposed_y_stand[indicator_imprecise, 3],
                   my_g = as.numeric(as.character(g)),
                   sd_of_y_precise = sd_of_y_precise,
                   rounding_degrees = rounding_degrees,
                   indicator_precise = indicator_precise,
                   indicator_imprecise = indicator_imprecise,
                   indicator_outliers = indicator_outliers,
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
      cat("Hessian matrix couldn't be inverted (in the imputation function of the rounded continuous variable).
              Still, you should get a result, but which needs special attention.\n")

      tmp <- diag(ncol(hess))
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
  if(ncol(X_model_matrix_1_sub) == 1){
    beta_hat <- matrix(1, ncol = 1)
  }else{
    beta_hat <- matrix(pars[grep("^coef_y_on_x", colnames(pars))], ncol = 1)
  }

  gamma1_hat <- pars[grep("^gamma1", colnames(pars))]
  if(length(gamma1_hat) == 0){
    gamma1_hat <- 0
  }
  gamma_hat <- matrix(pars[grep("^coef_g_on_psi", colnames(pars))], ncol = 1)
  sigma_hat <- pars[grep("^sigma", colnames(pars))]

  #For the potentially ROUNDED observations, the potential rounding tendency has to be modeled
  mu_g <- gamma1_hat * (as.matrix(Y_X_2_all[indicator_precise, xnames_1, drop = FALSE]) %*% beta_hat) +
    PSI_in_negloglik %*% gamma_hat

  #For the ROUNDED, MISSING and INTERVAL observations, y has to be modeled.
  mu_y <- as.matrix(Y_X_2_all[, xnames_1, drop = FALSE]) %*% beta_hat

  #so mu_g and mu_y generally differ in length!

  #The covariance matrix from equation (3)
  Sigma <- matrix(c(1 + gamma1_hat^2 * sigma_hat^2,
                      gamma1_hat * sigma_hat^2, gamma1_hat * sigma_hat^2,
                      sigma_hat^2), nrow = 2)

  ###########################################################
  #BEGIN IMPUTING INTERVAL-DATA AND COMPLETELY MISSING DATA#
  # The imputation for precise but rounded data follows in the next section.

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
  #do the backtransformation from
  #standardised to unstandardised
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
  actually_rounded <- g != 0
  y_lower <- (decomposed_y[indicator_precise, "precise"][actually_rounded] -
                half_interval_length[as.numeric(as.character(g))[actually_rounded]] - mean_of_y_precise)/sd_of_y_precise + 1

  y_upper <- (decomposed_y[indicator_precise, "precise"][actually_rounded] +
                half_interval_length[as.numeric(as.character(g))[actually_rounded]] - mean_of_y_precise)/sd_of_y_precise + 1

  g_upper <- bounds_for_g_hat[as.numeric(as.character(g))[actually_rounded] + 1]

  elements <- cbind(-Inf, mu_g[actually_rounded], g_upper,
                    y_lower, mu_y[indicator_precise, 1][actually_rounded], y_upper)

  # Note: we set g_lower to -Inf because we state that a value of 1500 is not necessarily
  # a multiple of 500; it could also be rounded to the next multiple of 10 or even 1.
  colnames(elements) <- c("g_lower", "mean_g","g_upper", "y_lower","mean_y",   "y_upper")

  ###indicator which of the precise observations need to be imputed due to rounding
  #(incl. rounding to the nearest number, i.e. rounding to degree 1)
  #(and not because they are missing)
  rounded <- rep(TRUE, sum(actually_rounded))

  counter <- 0
  while(any(rounded) & counter <= 1000){
    counter <- counter + 1

    if(counter == 1000){
      print(paste("Imputation of rounded continuous variable", colnames(y_df)[1], "failed."))
    }
    ###draw values for g and y from a truncated multivariate normal
    ###drawn y must be between y_lower and y_upper
    ###drawn g must be between g_lower and g_upper
    #If the model has some problems finding a suitable rounding degree,
    #the expected mean of g is increased to the threshold seperating the highest and second highest
    # category
    if(counter >= 100){
      elements[, 2] <- bounds_for_g_hat[length(bounds_for_g_hat)-1]
    }

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
    # We then check whether it is mu_g or mu_y, that lies outside its interval
    #and then replace the corresponding mean
    #by a uniform sample between the lower and the upper bound.


    #replace the invalid draws with valid ones.
    # For the latent rounding tendency, we use the highest possible rounding tendency
    # For y, we use a uniform sample between the highest and lowest possible
    #bounds of y.
    problematic_draws <- is.na(mytry[, 1])

    problematic_elements <- elements[problematic_draws, , drop = FALSE]

    #If the model has even further deficiancies, all remaining elements are considered problematic
    if(counter >= 200){
      toosmall_ys <- elements[, 5] < elements[, 4]
      toolarge_ys <- elements[, 5] > elements[, 6]
      elements[toosmall_ys, 5] <- elements[toosmall_ys, 4]
      elements[toolarge_ys, 5] <- elements[toolarge_ys, 6]
    }

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
    imp_tmp[indicator_precise][actually_rounded][rounded] <- imp_precise_temp

    #... but test if estimated rounding degree and proposed y can explain the observed y.
    # E.g. the estimated rounding degree 10 and the proposed y 2063 doesn't match
    #to an observed value 2100. A degree of 100 would match in this case.
    #If degree and y do match set the value for rounded to FALSE.
    # The remaining (non-matching) observations get a new proposal y and rounding degree.
    domatch <- floor(imp_precise_temp[, 1]/rounding_degrees[round_int] + 0.5) * rounding_degrees[round_int] ==
                          decomposed_y[indicator_precise, "precise"][actually_rounded][rounded]
    rounded[rounded][domatch] <- FALSE

  }

  #restore the original values for those observations with no rounding at all (rounding degree 0).
  no_rounding_at_all_clean <- no_rounding_at_all
  no_rounding_at_all_clean[is.na(no_rounding_at_all)] <- FALSE
  imp_tmp[no_rounding_at_all_clean] <- decomposed_y[no_rounding_at_all_clean, "precise"]
  y_ret <- data.frame(y_ret = imp_tmp)

  return(y_ret)
}
