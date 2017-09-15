

#' Sample imputation.
#'
#' Function to sample values in a variable from other (observed) values in this variable.
#' So this imputation does not use further covariates.
#' @param variable A vector of size \code{n} with missing values.
#' @return A n times 1 data.frame without missing values.
#' @examples
#' set.seed(123)
#' sample_imp(c(1, NA, 3, NA, 5))
#' @export
sample_imp <- function(variable){

  if(is.data.frame(variable)){
    stop("You passed a data.frame instead of a vector to sample_imp.")
  }
  if(all(is.na(variable))) stop("Variable consists only of NAs.")

  ret <- data.frame(target = variable)
  need_replacement <- is.na(variable) | is.infinite(variable)
  ret[need_replacement, 1] <- sample(size = sum(need_replacement),
                             variable[!need_replacement], replace = TRUE)
  return(ret)
}


#' Get the type of variables.
#'
#' Function checks whether a variable is: ...
#' \itemize{
#'  \item continuous (numeric values),
#'  \item semicontinuous (numeric values with more than 5\% of them are 0),
#'  \item rounded continuous (if continuous values are rounded to the closest multiple of 5, 10, 50, 100, 500 or 1000.
#'  We see this to be the case if more than 50\% of the observations are divisible by 5)
#'  \item count data (if all values are integers).
#'  \item an intercept (the same value for all observations),
#'  \item binary (two different values - like 0s and 1s or "m" and "f"),
#'  \item categorical (the variable is a factor or has more than 3 different values)
#'  \item ordered categorical (the categorical variable is ordered.)
#'}
#'
#' @param variable A variable (vector) from your data set.
#' @return A character denoting the type of \code{variable}.
#' @examples get_type(iris$Sepal.Length); get_type(iris$Species)
#' @export
get_type <- function(variable){

  if(all(is.na(variable))) return(NA)

  if(length(table(variable)) == 1){
    type <- "intercept"
    return(type)
  }

  if(length(table(variable)) == 2){
    type <- "binary"

    return(type)
  }


  if(is.integer(variable)){
    type <- "count"

    return(type)
  }

  if(is_interval(variable)){
    variable <- as.character(variable)
    type <- "interval"
    #if the precise part of an interval variable is rounded,
    #the whole variable is considered to be a rounded continous variable

    if(sum(decompose_interval(interval = variable)$precise %% 5 == 0, na.rm = TRUE)/
       length(variable[!is.na(variable)]) > 0.5){
      type <- "roundedcont"
    }
    return(type)
  }

  if(is.vector(variable) || is.factor(variable) || is.array(variable)){
    if(is.numeric(variable)){
      type <- "cont"

      # if there is no difference between a variable and their as.integer()-result,
      # it is considered to be an integer...
      if(max(abs(variable - as.integer(variable)), na.rm = TRUE) == 0){
        type <- "count"
        #... despite it is a rounded continuous variable:
        if(sum(variable %% 5 == 0, na.rm = TRUE)/
           length(variable[!is.na(variable)]) > 0.5){
          type <- "roundedcont"
        }
        return(type)
      }

      #if the variable is numeric and more than 5% of the variable values are 0,
      #We consider this variable as semi-continious.

      if(sum(variable == 0, na.rm = TRUE)/
         length(variable[!is.na(variable)]) > 0.05){
        type <- "semicont"
      }

      # if more than 50% of the data are divisible by 5, they are considered to be rounded
      # continous
      # Observed 0s shall not be counted as being divisible by 5,
      #because otherwise semi-continous variables might be considered to be rounded-continous.

      if((sum(variable %% 5 == 0, na.rm = TRUE) - sum(variable == 0, na.rm = TRUE))/
        length(variable[!is.na(variable)]) > 0.5){
        type <- "roundedcont"
      }

      return(type)

    }

    #checking the length of the table to be larger than two is not sufficient
    #to get a real categoriacal variable, because for a (real) continious variable
    #the length of the table is equal to the length of the variable
    #(what will be in general larger than two...)
    #Therefore it will be checked if the table has more than two entries AND
    #wheter the varriable is a factor.
    if(is.factor(variable) || length(table(variable)) > 2){
      type <- "categorical"
      if(is.ordered(variable)){
        type <- "ordered_categorical"
      }

      return(type)
    }

  }else{

    if(ncol(variable) == 1){
      ret <- get_type(variable[, 1])
      return(ret)
    }
    if(ncol(variable) > 1){
      if(is.matrix(variable)){
        variable <- as.data.frame(variable)
      }
      ret <- unlist(lapply(variable, get_type))

      return(ret)
    }
  }
}

#' Helps the user to make a list of types.
#'
#' In \code{hmi} the user can add a list of types. This function gives him a framework
#' with suggestions. Of course he can changes on his own afterwards.
#' For example, if a continuous variable as only two observations left, then get_type
#' interpret this as a binary variable and not a continuous.
#' @param data the data.frame also passed to \code{hmi}.
#' @return a list with suggested types. Each list element has the name of a variable
#' in the data.frame. The elements contain a single character denoting the type of the variable.
#' See \code{get_type} for details about the variable types.
#' @export
list_of_types_maker <- function(data){

  if(!is.data.frame(data)) stop("We need the data in the data.frame format!")
  ret <- lapply(data, get_type)

  return(ret)
}


#' Averages the results of the imputation function \code{hmi}.
#'
#' This function applies the analysis the user is interested in, on all different imputed dataset.
#' Then the results are pooled by simply averaging the results. So the user has to make sure that
#' his analysis produces results with a meaningful average. And furthermore has to accept that no
#' variance is calculated for these parameters.
#' @param mids A \code{mids} (multiply imputed data set) object.
#' Either from the \code{hmi} imputation function or \code{mice}.
#' @param analysis_function A user generated function that contains the model and the model parameters
#' he is interested in. See examples.
#' @return A vector with all averaged results.
#' @examples
#' \dontrun{
#' my.formula <- Reaction ~ Days + (1 + Days|Subject)
#' my_analysis <- function(complete_data){
#'  # In this list, you can write all the parameters you are interested in.
#'  # Those will be averaged.
#'  # So make sure that averaging makes sense and that you only put in single numeric values.
#'  parameters_of_interest <- list()
#'
#'  # ---- write in the following lines, what you are interested in to do with your complete_data
#'  # the following lines are an example where the analyst is interested in the fixed intercept
#'  # and fixed slope and the random intercepts variance,
#'  # the random slopes variance and their covariance
#'  my_model <- lmer(my.formula, data = complete_data)
#'
#'  parameters_of_interest[[1]] <- fixef(my_model)[1]
#'  parameters_of_interest[[2]] <- fixef(my_model)[2]
#'  parameters_of_interest[[3]] <- VarCorr(my_model)[[1]][1, 1]
#'  parameters_of_interest[[4]] <- VarCorr(my_model)[[1]][1, 2]
#'  parameters_of_interest[[5]] <- VarCorr(my_model)[[1]][2, 2]
#'  names(parameters_of_interest) <- c("beta_intercept", "beta_Days", "sigma0", "sigma01", "sigma1")
#'
#'  # ---- do not change this function below this line.
#'  return(parameters_of_interest)
#' }
#' require("lme4")
#' require("mice")
#' data(sleepstudy, package = "lme4")
#' test <- sleepstudy
#' test$Intercept <- 1
#' test[sample(1:nrow(test), size = 20), "Reaction"] <- NA
#' hmi_imp <- hmi(data = test, model_formula = my.formula, M = 5, maxit = 1)
#' hmi_pool(mids = hmi_imp, analysis_function = my_analysis)
#' #if you are interested in fixed effects only, consider using \code{pool} from \code{mice}.
#' pool(with(data = hmi_imp, expr = lmer(Reaction ~ Days + (1 + Days | Subject))))
#' }
#' @export
hmi_pool <- function(mids, analysis_function){

  if (!mice::is.mids(mids)){
    stop("Your multiple imputed data set must have class mids.")
  }

  results <- list()

  for (i in 1:mids$m) {
    results[[i]] <- analysis_function(mice::complete(mids, i))
  }

  tmp <- simplify2array(results)
  mode(tmp) <- "numeric"
  return(rowMeans(tmp))
}


#' calculate the likelihood contribution of the data
#'
#' This function based on Drechsler, Kiesl & Speidel (2015)
#' is needed in the imputation routine for rounded income.
#' It calculates the likelihood contribution of the data
#' (regardless whether they are observed precisely or presumably rounded).
#' @param para This is the vector \eqn{Psi} of parameters
#' (see p. 62 in Drechsler, Kiesl & Speidel, 2015).
#' With respect to them, the value returned by negloglik2 shall be
#' maximized.\cr
#' The starting values are c(kstart, betastart2, gammastart, sigmastart2)
#' (the 6 thresholds (or "cutting points") for the latent variable behind the rounding degree,
#' the regression parameters explaining the logged income,
#' the regression parameters explaining the rounding degree
#' and the variance parameter).
#' @param X the data.frame of covariates.
#' @param y_in_negloglik the target variable (a vector).
#' @param lower the lower bound of an interval variable.
#' @param upper the upper bound of an interval variable.
#' @param my_p This vector is the indicator of the (highest possible) rounding degree for an observation.
#' This parameter comes directly from the data.
#' @param mean_y_precise the scalar with the value of the mean of the target variable.
#' @param sd_y_precise the scalar with the value equal to the standard deviation of the target variable.
#' @references Joerg Drechsler, Hans Kiesl, Matthias Speidel (2015):
#' "MI Double Feature: Multiple Imputation to Address Nonresponse and Rounding Errors in Income Questions",
#' Austrian Journal of Statistics, Vol. 44, No. 2, \url{http://dx.doi.org/10.17713/ajs.v44i2.77}
#' @return An integer equal to the (sum of the) negative log-likelihood contributions (of the observations)
negloglik2 <- function(para, X, y_in_negloglik, lower = NA, upper = NA,
                       my_p, mean_y_precise, sd_y_precise){

  lower <- ifelse(is.na(lower), -Inf, lower)
  upper <- ifelse(is.na(upper), Inf, upper)

  # the first 6 parameters are the thresholds for the rounding degree
  # below k0: no rounding, between k0 and k1: rounding degree 5 etc.
  # The other degrees are 10, 50, 100, 500 and 1000
  k0 <- para[1]
  k1 <- para[2]
  k2 <- para[3]
  #k3 <- para[4]
  #k4 <- para[5]
  #k5 <- para[6]

  # the regression coefficients beta defining mean2 - the expected value of log(Y)
  # (see eq (2) in Drechsler, Kiesl & Speidel, 2015).
  # They also appear in mean1 (the expected value of G).
  reg_coefs <- matrix(para[4:(length(para) - 2)], ncol = 1)

  gamma1 <- para[length(para) - 1] #the gamma1 from eq (3)

  sigmax <- para[length(para)] # the variance for log(Y) (see eq (3) in Drechsler, Kiesl & Speidel, 2015)

  if (k1 < k0) return(1e50)        # make sure that threshholds
  if (k2 < k1) return(1e50)        # are in increasing order
  #if (k3 < k2) return(1e50)        #
  #if (k4 < k3) return(1e50)        #
  #if (k5 < k4) return(1e50)        #
  if (para[length(para)] <= 0) return(1e50) # the residual variance has to be positive
  n <- nrow(X)

  # mean and covariance of the joint normal of x and g ,
  # following Rubin/Heitjan

  #the mean for G (see eq (2) in Drechsler, Kiesl & Speidel, 2015)
  mean1 <- gamma1 * (as.matrix(X) %*% reg_coefs)

  #the mean for log(Y) (see eq (2) in Drechsler, Kiesl & Speidel, 2015)
  mean2 <- as.matrix(X) %*% reg_coefs

  #the covariance matrix for log(Y) and G (see eq (3) in Drechsler, Kiesl & Speidel, 2015)
  sigma <- matrix(c(1 + gamma1^2 * sigmax^2,
                    gamma1 * sigmax^2,
                    gamma1 * sigmax^2,
                    sigmax^2), nrow = 2)

  rho <- max(min(stats::cov2cor(sigma)[1, 2], 1), -1)

  sigma1 <- sqrt(sigma[1, 1])
  sigma2 <- sqrt(sigma[2, 2])

  # index the missing values
  missind_negloglik <- is.na(y_in_negloglik)
  y_obs <- y_in_negloglik[!missind_negloglik]

  #get the likelihood contributions of the imprecise observations
  interval_obs <- !is.na(lower) & !is.infinite(lower)
  sum_likelihood_imprecise_obs <- sum(log(contributions4intervals(lower = lower[interval_obs],
                                                                   upper = upper[interval_obs],
                                                                   mymean = mean2[interval_obs],
                                                                   mysd = sigma2)))

  # a1 are a7 the "components" of the log-liklihood contributions
  #
  # If the income is not divisible by 5, only a1 is relevant
  # if the income is divisible by 5, but not by 10, a1 + a2 are relevant
  # etc.

  #get the value of the standard bivariate normal cumulative density function
  #the better k0, k1, ..., sigma and rho are found, the better a1, a2, ...

  a1 <- pbivnormX(x =  c(k0 - mean1[!missind_negloglik])/sigma1,
                  y =  (y_obs + 0.5/sd_y_precise - mean2[!missind_negloglik])/sigma2,
                  rho = rho
  ) -
    pbivnormX(x =  c(k0 - mean1[!missind_negloglik])/sigma1,
              y =  (y_obs - 0.5/sd_y_precise - mean2[!missind_negloglik])/sigma2,
              rho = rho
    )

  a1 <- pmax(1e-100, abs(a1))
  # reason: if a1 is close to 0, a1 can be 0, due to computational imprecision


  a2 <- components(ki = k0, kj = k1,
                   mean1_obs = mean1[!missind_negloglik],
                   mean2_obs = mean2[!missind_negloglik],
                   sigma1 = sigma1,
                   sigma2 = sigma2,
                   rho = rho,
                   y_obs = y_obs,
                   mean_y_precise = mean_y_precise,
                   sd_y_precise = sd_y_precise,
                   half_divisor = 5)

  a3 <- components(ki = k1, kj = k2,
                   mean1_obs = mean1[!missind_negloglik],
                   mean2_obs = mean2[!missind_negloglik],
                   sigma1 = sigma1,
                   sigma2 = sigma2,
                   rho = rho,
                   y_obs = y_obs,
                   mean_y_precise = mean_y_precise,
                   sd_y_precise = sd_y_precise,
                   half_divisor = 50)

  a4 <- stats::pnorm(y_obs + 500/sd_y_precise, mean2[!missind_negloglik], sigma2) -

    pbivnormX( x =  c(k2 - mean1[!missind_negloglik])/sigma1,
               y =  (y_obs + 500/sd_y_precise - mean2[!missind_negloglik])/sigma2,
               rho = rho
    ) -

    stats::pnorm(y_obs - 500/sd_y_precise, mean2[!missind_negloglik], sigma2)  +

    pbivnormX( x =  c(k2 - mean1[!missind_negloglik])/sigma1,
               y =  (y_obs - 500/sd_y_precise - mean2[!missind_negloglik])/sigma2,
               rho = rho
    )

  a2 <- pmax(0, a2) # replacing negative values by 0.
  a3 <- pmax(0, a3)
  a4 <- pmax(0, a4)

  result_try <- cbind(a1 ,
                      a2 * (my_p[!missind_negloglik] >= 1),
                      a3 * (my_p[!missind_negloglik] >= 2),
                      a4 * (my_p[!missind_negloglik] >= 3))

  # Sum up the likelihood contributions (see eq (5) in Drechsler, Kiesl & Speidel, 2015)
  result <- log(rowSums(result_try, na.rm = TRUE))

  ret <- sum(result) + sum_likelihood_imprecise_obs
  return(-ret) # notice the minus
}

#' calculate the likelihood contribution of interval data only
#'
#' calculate the likelihood contribution of interval data only
#' @param para This is the vector \eqn{Psi} of parameters defining model
#' (see p. 62 in Drechsler, Kiesl & Speidel, 2015).
#' With respect to them, the value returned by this function shall be
#' maximized.\cr
#' The starting values are c(betastart2, sigmastart2)
#' (the regression parameters explaining the logged income,
#' and the variance parameter).
#' @param X the data.frame of covariates.
#' @param lower the lower bound of an interval variable.
#' @param upper the upper bound of an interval variable.
#' @references Joerg Drechsler, Hans Kiesl, Matthias Speidel (2015):
#' "MI Double Feature: Multiple Imputation to Address Nonresponse and Rounding Errors in Income Questions",
#' Austrian Journal of Statistics, Vol. 44, No. 2, \url{http://dx.doi.org/10.17713/ajs.v44i2.77}
#' @return An integer equal to the (sum of the) negative log-likelihood contributions (of the observations)
negloglik2_intervalsonly <- function(para, X, lower, upper){

  # the regression coefficients beta defining mean2 - the expected value of log(Y)
  # (see eq (2) in Drechsler, Kiesl & Speidel, 2015).
  reg_coefs <- matrix(para[1:(length(para) - 1)], ncol = 1)

  sigmax <- para[length(para)]
  # the variance for log(Y) (see eq (3) in Drechsler, Kiesl & Speidel, 2015)

  n <- nrow(X)
  result <- rep(NA, n)

  # mean and covariance of the joint normal of x and g ,
  # following Rubin/Heitjan

  #the mean for log(Y) (see eq (2) in Drechsler, Kiesl & Speidel, 2015)
  mean2 <- as.matrix(X) %*% reg_coefs

  #the covariance matrix for log(Y) and G (see eq (3) in Drechsler, Kiesl & Speidel, 2015)

  sigma2 <- sigmax

  #get the likelihood contributions of the imprecise observations
  interval_obs <- !is.na(lower)
  sum_likelihood_imprecise_obs <- sum(log(contributions4intervals(lower = lower[interval_obs],
                                                                   upper = upper[interval_obs],
                                                                   mymean = mean2[interval_obs],
                                                                   mysd = sigma2)))
  return(-sum_likelihood_imprecise_obs) # notice the minus
}

#' get the likelihood contributions of interval data
#'
#' This function calculates the likelihood contributions of interval data
#' @param lower a vector with the lower bounds of an interval covariate.
#' @param upper a vector with the upper bounds of an interval covariate.
#' @param mymean a numeric for the expected value of the normal distribution
#' (which is one of the parameters trying to be optimized
#' so that the likelihood becomes maximized)
#' @param mysd a numeric for the standard deviation of the normal distribution
#' (which is one of the parameters trying to be optimized
#' so that the likelihood becomes maximized)
#' @return a vector giving the likelihood contributions of the interval data.
contributions4intervals <- function(lower, upper, mymean, mysd){

  ret <- stats::pnorm(upper, mean = mymean, sd = mysd) -
    stats::pnorm(lower, mean = mymean, sd = mysd)

  return(ret)
}

#' calculate probabilities from the cumulative distribution function of a standard bivariate normal distribution
#'
#' A modified version of pbivnorm() from package \code{pbivnorm}.
#' It is needed in the imputation routine for rounded income.
#' @param x the vector (or a two columned matrix) with the values of the first random variable
#' @param y the vector with the values of the second random variable
#' @param rho the correlation (a scalar) between the two random variables.
#' @return A vector with the values of the density distribution at the points (x, y).
pbivnormX <- function (x, y, rho = 0) {

  #In case x is a matrix, take the second column as y and the first as x.
  if  (is.matrix(x)) {
    if (ncol(x) != 2)
      stop("'x' must have two columns if specified as a matrix")
    if (!missing(y) && !is.null(y)) #
      warning("'x' was specified as a matrix, so 'y' will be ignored")
    y <- x[, 2]
    x <- x[, 1]
  }

  if (any(abs(rho) > 1)) stop("'rho' must be a valid correlation (-1 <= rho <= 1)")

  if (length(x) != length(y)) stop("'x' and 'y' must have same length")

  x <- replace(x, x == Inf, 1e200)      #  function adjusted here
  x <- replace(x, x == -Inf, -1e200)    #  function adjusted here
  y <- replace(y, y == Inf, 1e200)      #  function adjusted here
  y <- replace(y, y == -Inf, -1e200)    #  function adjusted here

  ret <- pbivnorm::pbivnorm(x = x, y = y, rho = rho, recycle = TRUE)

  return(ret)
}


#' Function to get the likelihood contribution of different rounding degrees
#'
#' It is needed in the imputation routine for rounded income.
#' See equation (5) in Drechsler, Kiesl & Speidel (2015)
#' @param ki An integer for the i-th threshold (or "cutpoint")
#' @param kj An integer for the j-th threshold (or "cutpoint") (ki < kj)
#' @param mean1_obs A vector with the expected value of G for the observed data
#' @param mean2_obs A vector with the expected value of log(Y) for the observed data
#' @param sigma1 A scalar for the variance of the G
#' @param sigma2 A scalar for the variance of log(Y)
#' @param rho A scalar from [-1, 1] specifying the correlation between G and log(Y)
#' @param y_obs The vector of the target variable (with all NAs removed)
#' @param mean_y_precise A scalar specifying the mean of the target variable
#' @param sd_y_precise A scalar specifying the standard deviation of the target variable
#' @param half_divisor A scalar needed to find the bounds of possible rounding.
#' E.g. if rounding happens to the closest multiple of 500, the half.divisor is 250.
#' @references Joerg Drechsler, Hans Kiesl, Matthias Speidel (2015):
#' "MI Double Feature: Multiple Imputation to Address Nonresponse and Rounding Errors in Income Questions",
#' Austrian Journal of Statistics, Vol. 44, No. 2, \url{http://dx.doi.org/10.17713/ajs.v44i2.77}
#' @return A vector with the contribution to the likelihood of every individual with an observed target variable value
components <- function(ki, kj, mean1_obs, mean2_obs, sigma1, sigma2, rho, y_obs,
                       mean_y_precise, sd_y_precise, half_divisor){

  upper_inner <- ((y_obs + half_divisor - mean_y_precise)/sd_y_precise - mean2_obs)/sigma2
  lower_inner <- ((y_obs - half_divisor - mean_y_precise)/sd_y_precise - mean2_obs)/sigma2
  lower_outer <- c(ki - mean1_obs)/sigma1
  upper_outer <- c(kj - mean1_obs)/sigma1

  ret <- doubleintegral(lower_inner = lower_inner, upper_inner = upper_inner,
                        lower_outer = lower_outer, upper_outer = upper_outer,
                        cdf = pbivnormX, rho = rho)

  return(ret)
}

#' Function to calculate double integrals
#'
#' This function is primarily build to make the function \code{components} neater.
#' @param lower_inner The vector for the lower bound for the inner integral
#' @param upper_inner The vector for the upper bound for the inner integral
#' @param lower_outer The vector for the lower bound for the outer integral
#' @param upper_outer The vector for the upper bound for the outer integral
#' @param cdf the cumulative density function (from the class "function")
#' @param ... Further arguments passed to the cdf.
#' @return a vector with the value of the double integral for each observation (with an observed target variable)
doubleintegral <- function(lower_inner, upper_inner, lower_outer, upper_outer,
  cdf, ...){
  ret <- (
    cdf(x = upper_outer, y =  upper_inner, ...) -
      cdf(x = upper_outer, y =  lower_inner, ...)
  ) - (
    cdf(x = lower_outer, y =  upper_inner, ...)-
      cdf(x = lower_outer, y =  lower_inner, ...)
  )
  return(ret)
}

#' Function need to multivariate samples of a truncated multivariate normal distribution
#'
#' As rtmvnorm only allows one mean vector of one multivariate normal distribution,
#' but we need different mean vectors for different multivariate normal distributions,
#' we implement this function. This function in combination with \code{apply},
#' allows us to sample from a truncated multivariate normal distribution
#' with different mean vectors.
#' @param elements Originally a matrix, but when passed to samp, it is a vector.
#' The first length_mean elements are the mean vector of g and y,
#' the next two elements are the lower bounds for g and y,
#' the last two elements are the upper bounds for g and y.
#' @param Sigma The covariance matrix of the multivariate normal distribution to sample from.
#' @return A length_mean x 1 matrix with the samples for g and y.
sampler <- function(elements, Sigma){

  ret <- tmvtnorm::rtmvnorm(1,
                            mean = elements[1:2],
                            sigma = Sigma,
                            lower = elements[3:4],
                            upper = elements[5:6],
                            algorithm = "gibbs", burn.in.samples = 1000)
  return(ret)
}

#' Function to extract the different elements of a formula
#'
#' The function searches for the target variable, fixed effects variables,
#' if there is a cluster ID: this and the random effects variables.\cr
#' The names of the fixed and random intercepts variable (if existent) are explicitly labeled
#' In imputation models, the target variable can act as covariate
#' for other covariates - so we treat the target variable as fix effect variable.
#' @param model_formula A formula (from class \code{formula})
#' @param constant_variables A Boolean-vector of length equal to the number of columns in the data set
#' specifying whether a variable is a constant variable (eg. an intercept variable) or not.
#' @param variable_names_in_data A character-vector with the column names of the data set.
#' @return A list with the names of the target variable, the intercept variable,
#' the fixed and random effects covariates (which includes the name of the target variable),
#' the variables with interactions and the cluster id variable.\cr
#' If some of them don't exist, they get the value "".
#' @export
extract_varnames <- function(model_formula = NULL, constant_variables,
                             variable_names_in_data){


  # Set up default values for key variables like the target variable, the clusterID,
  # and the random covariates
  target_varname <- NULL

  #extract key variables from formula
  if(!is.null(model_formula)){

    # if it was a character string, make it a formula
    if(class(model_formula) == "character"){
      warning("We need your model_formula to be a formula (and not a character).
              So we changed it automatically, but it would be better if you do it.")

      model_formula <- stats::formula(model_formula)
    }

    variable_names_in_formula <- all.vars(stats::terms(model_formula))

    #looking for possible interactions
    interactions <- grep(":", x = attr(stats::terms(model_formula), "term.labels"),
                         value = TRUE)
    if(length(interactions) == 0){
      interactions <- ""
    }

    # ----------Target variable
    target_varname_full <- as.character(model_formula)[2]

    #check if the formula contains any functions like "log(y)" in log(y) ~ x.
    target_varname <- gsub(".*\\((.*)\\).*", "\\1", target_varname_full)

    if(target_varname != target_varname_full){
      warning(paste("For the tool to work saver we suggest to do the transformation -->", target_varname_full,
                    "<-- before running the imputation."))
    }

    # -------- Cluster ID

    clID_varname <- sapply(lme4::findbars(model_formula), function(x) as.character(x)[3])

    #check if there is a cluster variable specified:...
    if(any(is.na(clID_varname)) | length(clID_varname) == 0){
      #... if it is not, we don't have a random effects model
      clID_varname <- ""
      random_intercept_exists <- FALSE
      randomeffects_varname <- ""

    }else{#if the cluster id is not NA


      # ----------- Z
      nearly_randomeffects_varname <-
        strsplit(as.character(lme4::findbars(model_formula)), "\\|")[[1]][1]

      # split the variables up
      nearly_randomeffects_varname <-
        strsplit(nearly_randomeffects_varname, "(\\+|\\*|\\:)")[[1]]

      # remove spaces
      randomeffects_varname <- gsub(" ", "", nearly_randomeffects_varname)

      # -- check for random intercept
      # If not explicitely removed a random intercept is always present
      # cf. lmer(mpg ~ cyl + (drat|am), data = mtcars)
      random_intercept_exists <- TRUE
      # Only if a "0" or "-1" is in the model_formula, no random intercept is estimated
      # cf. lmer(mpg ~ cyl + (0 + drat|am), data = mtcars)

      if("0" %in% randomeffects_varname | "-1" %in% randomeffects_varname){
        random_intercept_exists <- FALSE

        #exclude random effects called "0"
        randomeffects_varname <- randomeffects_varname[randomeffects_varname != "0"]

        #exclude random effects called "-1"
        randomeffects_varname <- randomeffects_varname[randomeffects_varname != "-1"]
      }

      # To be fool proof: if the user puts a 1 in his model_formula,
      # a random intercept shall be calculated, even if the model could do something different
      # e.g. in the case lmer(mpg ~ cyl + (1 + 0 + drat|am), data = mtcars)
      if("1" %in% randomeffects_varname) random_intercept_exists <- TRUE

    }

    # ---------X variables

    fixedeffects_varname <- all.vars(stats::delete.response(stats::terms(lme4::nobars(model_formula))))
    # --check if intercept is present
    fixed_intercept_exists <- attributes(stats::delete.response(stats::terms(lme4::nobars(model_formula))))$intercept == 1

    # --remove cluster ID from the X-variables
    # (but, if model_formula was specified correctly it shouldn't be there anyway)
    fixedeffects_varname <- fixedeffects_varname[ ! fixedeffects_varname %in% clID_varname]


    if("0" %in% fixedeffects_varname | "-1" %in% fixedeffects_varname){
      fixed_intercept_exists <- FALSE

      #exclude fixed effects called "0"
      fixedeffects_varname <- fixedeffects_varname[fixedeffects_varname != "0"]

      #exclude fixed effects called "-1"
      fixedeffects_varname <- fixedeffects_varname[fixedeffects_varname != "-1"]
    }

    ## Handling of a intercept variable.

    # If the model_formula includes a fix or random intercept
    # or if there is one constant variable in the data set,
    # get a good name for this intercept variable and make sure that an
    # intercept variable with this name will exist in the data set.
    intercept_varname <- ""
    if(fixed_intercept_exists | random_intercept_exists | sum(constant_variables) == 1){

      # give intercept_varname a default value.
      intercept_varname <- "Intercept"

      # but if there is a constant variable, its name will be used  the intercept variable name
      if(sum(constant_variables) == 1){
        intercept_varname <- variable_names_in_data[constant_variables]
      }


      if(fixed_intercept_exists){

        # Make the fixed intercept a part of the fixed effects.
        fixedeffects_varname <- c(intercept_varname, fixedeffects_varname)

        # If the intercept_varname is not "1"...
        if(intercept_varname != "1"){
          #... then every fixedeffects_varname that is "1" has to be removed
          fixedeffects_varname <- fixedeffects_varname[ !fixedeffects_varname == "1"]
        }
      }


      # If the intercept_varname is not "1",
      # then every randomeffects_varname that is "1" has to be removed
      if(random_intercept_exists){

        randomeffects_varname <- c(intercept_varname, randomeffects_varname)

        if(intercept_varname != "1"){
          randomeffects_varname <- randomeffects_varname[ !randomeffects_varname == "1"]
        }
      }


      # Replace the "1" in fixedeffects_varname by the name of the fixed intercept variable.
      fixedeffects_varname[fixedeffects_varname == "1"] <- intercept_varname

      # remove doublets
      fixedeffects_varname <- unique(fixedeffects_varname)

      # Replace the "1" in randomeffects_varname by the name of the intercept variable.
      randomeffects_varname[randomeffects_varname == "1"] <- intercept_varname
      randomeffects_varname <- unique(randomeffects_varname)

    }


  }else{ # if model_formula is NULL
    # in this case impute every variable based on the others in a single level framework

    intercept_varname <- names(constant_variables)[constant_variables]
    fixedeffects_varname <- variable_names_in_data
    randomeffects_varname <- ""
    interactions <- ""
    clID_varname <- ""
  }

  # Note: there can be only one intercept variable in the data set.

  ret <- list(target_varname = target_varname,
              intercept_varname = intercept_varname,
              fixedeffects_varname = c(target_varname, fixedeffects_varname),
              randomeffects_varname = randomeffects_varname,
              interactions = interactions,
              clID_varname = clID_varname)
  return(ret)
}

#' Function to check multilevel models on the existence of fixed intercepts
#'
#' Function to check multilevel models on the existence of fixed intercepts.
#' The specification of an intercept by calling a 1-column (eg "int")
#' is not counted towards the existence of an intercept.
#' Contradictionary inputs like "~ 1 + 0 + X1 + ..." or "~ -1 + 1 + X1 + ..."
#' will throw an error.
#' @param model_formula A formula (from class \code{formula})
#' @return A boolean value indicating whether there is a fixed intercept in the model or not
#' @export
fixed_intercept_check <- function(model_formula){

  if(!is.null(model_formula)){

    # if it was a character string, make it a formula
    if(class(model_formula) == "character"){
      warning("We need your model_formula to be a formula (and not a character).
              So we changed it automatically, but it would be better if you do it.")

      model_formula <- stats::formula(model_formula)
    }
  }

  fixedeffects_varname <- all.vars(stats::delete.response(stats::terms(lme4::nobars(model_formula))))
  # --check if intercept is present
  fixed_intercept_exists <- attributes(stats::delete.response(stats::terms(lme4::nobars(model_formula))))$intercept == 1

  #check on unplausibilities:
  if(fixed_intercept_exists &
     ("0" %in% fixedeffects_varname | "-1" %in% fixedeffects_varname)){
    stop("Your model_formula specify both, the existence and absence of a fixed intercept.")
  }
  return(fixed_intercept_exists)
}


#' Function to check multilevel models on the existence of random intercepts
#'
#' Function to check multilevel models on the existence of random intercepts.
#' The specification of an intercept by calling a 1-column (eg "int")
#' is not counted towards the existence of an intercept.
#' Contradictionary inputs like "~ 1 + 0 + X1 + ..." or "~ -1 + 1 + X1 + ..."
#' will throw an error.
#' @param model_formula A formula (from class \code{formula})
#' @return A boolean value indicating whether there is a fixed intercept in the model or not
#' @export
random_intercept_check <- function(model_formula){

  if(!is.null(model_formula)){

    # if it was a character string, make it a formula
    if(class(model_formula) == "character"){
      warning("We need your model_formula to be a formula (and not a character).
              So we changed it automatically, but it would be better if you do it.")

      model_formula <- stats::formula(model_formula)
    }
  }

  randomeffects_varname <- as.character(lme4::findbars(model_formula)[[1]][2])
  # --check if intercept is present
  #supported ways to specify an intercept:
  # "1 +..."
  # "1+..."
  # "...+1"
  # "...+ 1"
  # "..."
  pro_int <- length(grep("^1\\ +|^1\\+|\\+1|\\+ 1", randomeffects_varname)) > 0
    #supported ways to supress an intercept
  # "-1 +..."
  # "- 1 +.."
  # "-1+..."
  # "- 1+.."
  # "...- 1"
  # "...-1"
  # "0 +..."
  # "0+.."
  # "...+ 0"
  # "...+0"
  contra_int <-
    length(grep("\\-1\\ +|\\- 1\\ +|\\-1\\+|\\- 1\\+|\\- 1|\\-1|^0\\ +|^0\\+|\\+ 0|\\+0",
                  randomeffects_varname)) > 0


  random_intercept_exists <- !contra_int
  if(pro_int & contra_int){
    print("We found indications for and against a random intercept in your model_formula.")
  }

  return(random_intercept_exists)
}

#' Function to transform numeric (or character) vectors into an interval object
#'
#' Function to transform numeric (or character) vectors into an interval object
#' @param x An object to transform.
#' Currently the function can transform numeric vectors and characters
#' @seealso \link[hmi]{generate_interval}
#' @return A vector of class \code{interval}.
#' @examples as.interval(c("1000;2000", "700:700"))
#' as.interval(c("1500;1500", "700:700"))
#' @export
as.interval <- function(x){
  #Therefore it is needed, that "NA;NA" is a level in interval.
  if(is.factor(x)){
    x <- as.character(x)
  }
  if(is.numeric(x)){
    x <- paste(x, ";", x, sep = "")
  }
  x[is.na(x)] <- c("NA;NA")
  raw <- unlist(strsplit(x, ";"))
  raw[raw == "NA"] <- NA
  lower <- as.numeric(raw[seq(1, 2*length(x), by = 2)])
  upper <- as.numeric(raw[seq(2, 2*length(x), by = 2)])
  ret <- generate_interval(lower, upper)
  class(ret) <- "interval"
  return(ret)
}

#' Function to generate an interval object
#'
#' Function to generate transform into an
#' \code{interval} object from numeric (or character) vectors.
#' @param lower a vector with the lower values of a variable
#' @param upper a vector with the upper values of a variable
#' @return a character vector with the lower and upper bound values.
#' @export
generate_interval <- function(lower, upper){
  if(length(lower) != length(upper)) stop("lower and upper need the same length")

  # NA get the most uniformative values:
  # for the lower bound it is -Inf, and for the upper it is +Inf

  lower[is.na(lower)] <- -Inf
  upper[is.na(upper)] <- Inf


  lower[is.nan(lower)] <- -Inf
  upper[is.nan(upper)] <- Inf

  if(any(lower > upper)){
    stop("each element in lower must be smaller or equal to the corresponding element in upper")
  }
  ret <- paste(lower, upper, sep = ";")
  class(ret) <- "interval"
  return(ret)
}


#' Function to check whether an object is an interval
#'
#' If there are numerics separated by semicolons (;), this is considered to be an interval.
#' intervals with 2.4e5 are not considered to be an interval.
#' @param x an object
#' @return a boolean value indicting whether x is an interval or not
is_interval <- function(x){
  if("interval" %in% class(x)) return(TRUE)

  #NAs in x shall be ignored.
  x <- x[!is.na(x)]

  #search at the start for one or more digits optionally preceded by an "-"
  # optionally followed by one period (or "full stop") and one or more digits
  #or the interval starts with "-Inf"
  reg_low <- "^((-){0,1}[[:digit:]]{1,}(\\.[[:digit:]]{1,}){0,1}|-Inf|NA)"

  #the end of the interval is one or more digits optionally preceded by an "-"
  # optionally followed by one period (or "full stop") and one or more digits
  #or the interval ends with "Inf"
  reg_up <- "((-){0,1}[[:digit:]]{1,}(\\.[[:digit:]]{1,}){0,1}|Inf|NA)$"
  matches <- grep(paste(reg_low, reg_up, sep = ";"), x)

  return(length(matches) == length(x))

}

#' Transform interval objects into data.frames
#'
#' This function transforms interval objects into data.frames.
#' This is not only relevant on its own, it is also needed whenever
#' a function need objects as a data.frame (e.g. \code{View} or \code{cbind}).
#' @param x an object of class interval.
#' @param ... arguments to be passed to \code{as.data.frame}.
#' @return a data.frame containing x as a character
as.data.frame.interval <- function(x, ...){
  tmp <- as.data.frame(as.character(x), ...)
  colnames(tmp) <- colnames(x)
  return(tmp)
}

#' Split up intervals
#'
#' This function splits an interval object up into the lower and upper bound
#' @param interval an interval object of length n
#' (if it is something else, it is returned unchanged)
#' @return a n times 2 matrix. The first column is the lower bound, the second the upper bound.
#' @export
split_interval <- function(interval){

  #if it is no interval object, just return the object
  if(!is_interval(interval)) return(interval)

  #NA values have to be set temporarily to "NA;NA"
  #Therefore it is needed, that "NA;NA" is a level in interval.
  if(is.factor(interval)){
    interval <- as.character(interval)
  }
  interval[is.na(interval)] <- c("NA;NA")
  raw <- unlist(strsplit(interval, ";"))
  raw[raw == "NA"] <- NA
  lower <- as.numeric(raw[seq(1, 2*length(interval), by = 2)])
  upper <- as.numeric(raw[seq(2, 2*length(interval), by = 2)])
  return(matrix(c(lower, upper), ncol = 2))
}

#' decompose up intervals
#'
#' This function decomposes an interval object up into precise observations
#' (e.g. "1850.23;1850.23" into 1850.23),
#' imprecise observations (e.g. "1800;1900") and
#' missing observations ("-Inf;Inf" into NA)
#' @param interval an \code{interval} object of length n
#' (if it is something else, it is returned unchanged)
#' @return A data.frame with 5 columns.
#' 1. A column "precise" for the precise observations
#' (length of interval = 0, e.g. "3000;3000"). If observation i is not precise,
#' the i-th entry in this columns will be \code{NA}.
#' c("2500;2500", "3000;4000", "-Inf;0", NA) will lead to c(2500, NA, NA, NA)
#' 2. A column "lower" for the values of the lower bounds
#' of the imprecise observations (length of interval > 0,
#' e.g. "3000;4000" or "-Inf;0"), precise observations will get
#' \code{NA}s here.
#' c("2500;2500", "3000;4000", "-Inf;0", NA) will lead to c(NA, 3000, -Inf, NA)
#' 3. A column "upper" for the values of the upper bounds
#' of the imprecise observations.
#' c("2500;2500", "3000;4000", "-Inf;0", NA) will lead to c(NA, 4000, 0, NA)
#' 4. A column "lower_general" with the lower bound values of all observations,
#' without distinction between precise, imprecise or missing observations.
#' c("2500;2500", "3000;4000", "-Inf;0", NA) will lead to c(2500, 3000, -Inf, -Inf)
#' 5. A column "upper_general" with the upper bound values of all observations.
#' c("2500;2500", "3000;4000", "-Inf;0", NA) will lead to c(2500, 4000, 0, Inf)
decompose_interval <- function(interval){
  if(is.factor(interval)){
    interval <- as.character(interval)
  }
  if(!is_interval(interval)){
    ret <- data.frame(precise = interval, lower_imprecise = NA, upper_imprecise = NA,
                     lower_general = ifelse(is.na(interval), -Inf, interval),
                     upper_general = ifelse(is.na(interval), Inf, interval))
    return(ret)
  }
  tmp <- split_interval(interval)
  dists <- tmp[, 2] - tmp[, 1]
  precise <- ifelse(abs(dists) < 1e-20, tmp[, 1], NA)
  imprecise <- ifelse(abs(dists) >= 1e-20 | is.na(dists), interval, NA)
  tmp_imprecise <- split_interval(imprecise)
  ret <- data.frame(precise, lower_imprecise = tmp_imprecise[, 1],
                    upper_imprecise = tmp_imprecise[, 2],
                    lower_general = ifelse(is.na(tmp[, 1]), -Inf, tmp[, 1]),
                    upper_general = ifelse(is.na(tmp[, 2]), Inf, tmp[, 2]))
  return(ret)
}

#' Get standard NAs from interval data
#'
#' This function replaces observations with "-Inf;Inf" with the standard NAs (therefore 'sna')
#' @param x can by any object, but the function was designed for \code{interval}-objects.
#' @return In case of \code{x} being an \code{interval}-object, it returns a n times 2 matrix.
#'  The first column is the lower bound, the second the upper bound.
#'  Otherwise it returns just \code{x}.
#' @export
sna_interval <- function(x){

  #if it is no interval object, just return the object
  if(!is_interval(x)){
    return(x)
  }

  #replace observations with "-Inf;Inf" by NA
  x[x == "-Inf;Inf"] <- NA
  return(x)
}

#' Transform interval variables to an interval data frame
#'
#' This function is the path from this hmi package to the linLIR package (Wiencierz, 2012).
#' @param interval an \code{interval}
#' @return an interval data frame (idf-object) with one variable (having a lower and an upper bound).
#' @export
interval2idf <- function(interval){
  #lapply(interval, split_interval)
  dat <- as.data.frame(split_interval(interval))
  ret <- linLIR::idf.create(dat)
  return(ret)
}


#' Transform interval data frames into data.frames with interval variables
#'
#' This function is the path from the linLIR package (Wiencierz, 2012) to this hmi package.
#' @param idf an interval data frame (idf-object).
#' @return interval an \code{interval}.
#' @export
idf2interval <- function(idf){
  ret <- data.frame(matrix(nrow = idf$n, ncol = length(idf) - 1))
  for(i in 1:ncol(ret)){
    ret[, i] <- generate_interval(lower = idf[[i]][, 1], upper = idf[[i]][, 2])
  }
  colnames(ret) <- names(idf)[-length(idf)]
  return(ret)
}

#' Adding function
#'
#' Function to add single elements or vectors (of correct dimension) to the interval object
#' @param interval an object from class interval
#' @param x an single element or vector to add to the interval object
#' @return an interval object
#' @export
#' @rdname interval-add
"+.interval" <- function(interval, x){
  tmp <- split_interval(interval) + split_interval(x)
  generate_interval(lower = tmp[, 1], upper = tmp[, 2])
}


#' Subtraction function
#'
#' Function to subtract single elements or vectors (of correct dimension) from the interval object
#' @param interval an object from class interval
#' @param x an single element or vector to subtract to the interval object
#' @return an interval object
#' @export
#' @rdname interval-subtract
"-.interval" <- function(interval, x){
  tmp <- split_interval(interval) - split_interval(x)
  generate_interval(lower = tmp[, 1], upper = tmp[, 2])
}


#' Multiplication function
#'
#' Function to multiply single elements or vectors (of correct dimension) to the interval object
#' @param interval an object from class interval
#' @param x an single element or vector for multiplication
#' @return an interval object
#' @export
#' @rdname interval-multiply
"*.interval" <- function(interval, x){
  tmp <- split_interval(interval) * split_interval(x)
  generate_interval(lower = tmp[, 1], upper = tmp[, 2])
}

#' Dividing function
#'
#' Function to divide single elements or vectors (of correct dimension) to the interval object
#' @param interval an object from class interval
#' @param x an single element or vector for division
#' @return an interval object
#' @export
#' @rdname interval-divide
"/.interval" <- function(interval, x){
  tmp <- split_interval(interval) / split_interval(x)
  generate_interval(lower = tmp[, 1], upper = tmp[, 2])
}

#' Function to give the center of the interval
#'
#'  Function to give the center of the interval object
#' @param interval an object from class interval
#' @param inf2NA logical. If \code{TRUE}, entries containing -Inf or Inf, will return NA.
#' @return A numeric vector
#' @export
center.interval <- function(interval, inf2NA = FALSE){
  if(!is_interval(interval)){
    return(interval)
  }
  tmp <- rowMeans(x = split_interval(interval))
  if(inf2NA){
    tmp[is.infinite(tmp)] <- NA
  }
  return(tmp)
}

#' is.na for interval objects
#'
#' This functions checks whether elements from an \code{interval} object are \code{NA}
#' @param interval An \code{interval} object of length n
#' @return A boolean vector of length n indicating whether the entries in \code{interval}
#' are \code{NA} or not. Cf. \code{is.na}.
is.na.interval <- function(interval){
  ret <- is.na(as.character(sna_interval(interval)))
  return(ret)
}


#' Standardizing function
#'
#' Function to standardize variables that are numeric (continuous and count variables) but no rounded continuous, semicontinuous, intercepts or categorical variables.
#' @param X A n times p data.frame with p fixed (or random) effects variables.
#' @return A n times p data.frame with the standardized versions of the numeric variables.
#' @export
stand <- function(X){
  if(!is.data.frame(X)) stop("X has to be a data.frame.")
  types <- array(NA, dim = ncol(X))
  for(i in 1:length(types)){
    types[i] <- get_type(X[, i, drop = FALSE])
  }
  need_stand_X <- types %in% c("cont", "count")
  X_stand <- X
  X_stand[, need_stand_X] <- scale(X[, need_stand_X])
  return(X_stand)
}

#' Removing excessive factors
#'
#' Function to exclude variables that have too many levels as they may cause numerical problems
#' @param X A n times p data.frame with p fixed (or random) effects variables.
#' @param k An integer defining the allowed maximum of levels in a factor covariate.
#' @return A n times (p-r) data.frame, with r being the number of variables with too many factors.
#' @export
remove_excessives <- function(X, k = 10){

  if(!is.data.frame(X)) stop("X has to be a data.frame.")

  if(!is.numeric(k)) stop("k has to be a single numeric value.")
  if(length(k) > 1)  stop("k has to be a single numeric value.")

  # work out the type of the variables
  types <- array(dim = ncol(X))
  for(j in 1:length(types)) types[j] <- get_type(X[, j])

  # currently, R interprets interval data as categorical.
  # This can be unfeasible as there could be as many categories as observations.
  # Two practical approaches are 1. ignoring interval data. 2. splitting them up into
  # two separate variables. Currently we chose option 2.
  interval <- types == "interval"
  added_variables <- as.data.frame(matrix(NA, nrow = nrow(X), ncol = 2 * sum(interval)))
  interval_counter <- 1
  for(l in which(interval)){
    lower <- sample_imp(split_interval(X[, l])[, 1])
    upper <- sample_imp(split_interval(X[, l])[, 2])
    added_variables[, 2 * (interval_counter - 1) + 1:2] <- cbind(lower, upper)

    interval_counter <- interval_counter + 1
  }

  # determine the "categorical" variables.
  categorical <- types == "categorical"

  # get the name of the variables, that have too many levels
  too_many_levels <- colnames(X[, categorical, drop = FALSE])[
    apply(X[, categorical, drop = FALSE], 2, function(x) nlevels(factor(x))) > k]

  tmp <- X[, !names(X) %in% too_many_levels & !interval, drop = FALSE]

  ret <- cbind(tmp, added_variables)
  if(ncol(ret) == 0) stop("No variables left after removing excessive variables.")

  return(ret)
}


#' Cycling
#'
#' Function to do one imputation cycle on the given data. The function cycles through
#' every variable sequentially imputing the values, that are NA in the original data set
#' in that current variable. The function determines the type of the variable
#' and calls the suitable imputation function.
#' @param data_before The data.frame with the variables to impute.
#' @param original_data The original data.frame the user passed to \code{hmi}.
#' @param fe A list with the decomposed elements of the \code{model_formula}.
#' @param interaction_names A list with the names of the variables
#' that have been generated as interaction variables
#' @param NA_locator A n x p matrix localizing the missing values in the original
#' dataset. The elements are TRUE if the original data are missing and FALSE if the
#' are observed.
#' @param list_of_types a list specifying the types of the variables.
#' See \code{hmi} for details.
#' @param nitt An integer defining number of MCMC iterations (see \code{MCMCglmm}).
#' @param thin An integer defining the thinning interval (see \code{MCMCglmm}).
#' @param burnin An integer defining the percentage of draws from the gibbs sampler
#' that should be discarded as burn in (see \code{MCMCglmm}).
#' @param mn An integer defining the minimum number of individuals per cluster.
#' @return A data.frame where the values, that have a missing value in the original
#' dataset, are imputed.
#' @export
imputationcycle <- function(data_before,
                            original_data,
                            NA_locator,
                            fe, interaction_names,
                            list_of_types,
                            nitt, thin, burnin, mn){

  if(!is.data.frame(data_before)){
    stop("data_before has the be a data.frame")
  }

  if(!is.matrix(NA_locator)){
    stop("NA_locator has to be a matrix")
  }

  original_ordering <- colnames(data_before)
  missing_rates <- colMeans(NA_locator)

  #Sort data by missing rate
  data_before <- data_before[, names(sort(missing_rates)), drop = FALSE]
  NA_locator <- NA_locator[, names(sort(missing_rates)), drop = FALSE]

  #get variables with missing values
  incomplete_variables <- names(missing_rates)[missing_rates > 0]

  variables_to_impute <- union(incomplete_variables, names(list_of_types)[list_of_types == "roundedcont" |
                           list_of_types == "interval"])

  for(l2 in variables_to_impute){#do the imputation cycle

    #restore the original status of the variable which is to be imputed
    data_before[, l2] <- original_data[, l2]

    # get type
    tmp_type <- list_of_types[l2][[1]]

    #update the interaction variables
    if(fe$interactions != ""){
      for(j in 1:length(fe$interactions)){
        #get the product of the interaction variables
        interaction <- apply(data_before[, strsplit(fe$interactions[j], ":")[[1]], drop = FALSE], 1, prod)

        #add the interaction to the dataset
        data_before[, interaction_names[[j]]] <- interaction

        #This step has to be repeated after each imputation of a variable.
      }
    }else{
      interaction_names <- ""
    }

    #generate temporary data set, that contains only the variable to impute and
    #the covariates where no missing values occurs

    # First X and Z are generated
    # X contains the fixed effects variables and interaction variables
    X <- data_before[, union(fe$fixedeffects_varname[fe$fixedeffects_varname %in% names(data_before)],
                             interaction_names[[1]][interaction_names[[1]] %in% names(data_before)]),
                 drop = FALSE]

    # now safe the variables used for modeling random effects
    Z <- data_before[, fe$randomeffects_varname[fe$randomeffects_varname %in% names(data_before)],
                 drop = FALSE]

    tmp_X <- X[, names(X) != l2, drop = FALSE]
    tmp_Z <- Z[, names(Z) != l2, drop = FALSE]

    #now, remove variables from X an Z that currently have NAs
    #(wich should only occur in the first, initial imputation cycle)
    tmp_X <- tmp_X[, stats::complete.cases(t(tmp_X)), drop = FALSE]
    tmp_Z <- tmp_Z[, stats::complete.cases(t(tmp_Z)), drop = FALSE]

    #To preserve the relation of y and x properly,
    #y has to be included as covariate in the imputation model for x.
    #But not only as a fix effect covariate,
    #but also as a random effect covariate in the cases where x
    #is a random effects covariate in the analysis model for y.
    if(l2 %in% fe$randomeffects_varname){
      tmp_Z <- cbind(tmp_Z, X[, fe$target_varname])
    }

    #check number of observed values in each cluster.

    if(fe$clID_varname != ""){

      # Make sure that there are no missing values in the cluster ID:
      data_before[, fe$clID_varname] <- sample_imp(data_before[, fe$clID_varname])

      tab_1 <- table(data_before[!is.na(data_before[, l2]), fe$clID_varname])

      #merge small clusters to a large big cluster
      safety_counter <- 0
      while(any(tab_1 < mn) & safety_counter < 100){
        safety_counter <- safety_counter + 1

        #step 1: determine the smallest and the second smallest cluster.
        smallest_cluster <- tab_1[tab_1 == sort(as.numeric(tab_1))[1]][1]
        tab_1_without_smallest_cluster <- tab_1[names(tab_1) != names(smallest_cluster)]
        second_smallest_cluster <- tab_1_without_smallest_cluster[
                                   tab_1_without_smallest_cluster ==
                                  sort(as.numeric(tab_1_without_smallest_cluster))[1]][1]

        #step 2: new cluster name
        new_name <-  names(second_smallest_cluster)

        levels(data_before[, fe$clID_varname])[levels(data_before[, fe$clID_varname]) ==
                                          names(smallest_cluster)] <- new_name
        levels(data_before[, fe$clID_varname])[levels(data_before[, fe$clID_varname]) ==
                                          names(second_smallest_cluster)] <- new_name
        #step 3: repeat.

        tab_1 <- table(data_before[!is.na(data_before[, l2]), fe$clID_varname])

      }

    }


    # if too few observations are left for a usefull imputation,
    # the function stops
    if(sum(!is.na(data_before[, l2])) <= 2){
      stop(paste("Too few observations left in variable", l2))
    }

    #If no complete covariates are left
    #(which currently should not happen, as we currently force X to include an intercept variabel)
    #We have to do sample imputation.
    if(ncol(tmp_X) == 0){
      tmp_type <- "sample_imp"
      imp <- sample_imp(data_before[, l2])
    }

    #check whether a single level or multilevel imputation should be run
    # the standard is basically a multilevel imputation, but if one (or more) of the following
    # conditions is met, a single level impuation is run:
    # 1. there was no cluster ID found by the extract_varnames function
    # 2. the cluster ID cannot be found in the data.
    # 3. the variable to impute is the cluster ID.
    # 4. tmp_Z is currently empty
    # Explanation for 1 and 2: Without an cluster ID in the formula or the data, no multilevel model.
    # Explanation for 3: multilevel models say that the influence of the covariates changes from cluster to cluster
    # so cluster IDs are part of the covariates and not the target variable.
    # Explanation for 4: It can happen that the only random effects variable is going to be imputed,
    # resulting in an empty tmp_Z. But without random effects variables, no random effects model.
    use_single_level <- fe$clID_varname == "" || !(fe$clID_varname %in% names(data_before)) ||
      ncol(tmp_Z) == 0

    if(tmp_type == "binary"){
      if(use_single_level){

        #print("Currently the single level model.")
        imp <- imp_binary_single(y_imp = data_before[, l2],
                                 X_imp = tmp_X)

      }else{

        #print("Currently the multilevel model.")
        imp <- imp_binary_multi(y_imp = data_before[, l2],
                                X_imp = tmp_X,
                                Z_imp = tmp_Z,
                                clID = data_before[, fe$clID_varname],
                                nitt = nitt,
                                thin = thin,
                                burnin =  burnin)
      }
    }

    if(tmp_type == "cont"){

      if(use_single_level){

        #print("Currently the single level model.")
        imp <- imp_cont_single(y_imp = data_before[, l2],
                               X_imp = tmp_X)

      }else{

        #print("Currently the multilevel model.")
        imp <- imp_cont_multi(y_imp = data_before[, l2],
                              X_imp = tmp_X,
                              Z_imp = tmp_Z,
                              clID = data_before[, fe$clID_varname],
                              nitt = nitt,
                              thin = thin,
                              burnin = burnin)

      }

    }

    if(tmp_type == "semicont"){


      if(use_single_level){

        #print("Currently the single level model.")
        imp <- imp_semicont_single(y_imp = data_before[, l2],
                                   X_imp = tmp_X,
                                   heap = 0)
      }else{

        #print("Currently the multilevel model.")
        imp <- imp_semicont_multi(y_imp = data_before[, l2],
                                  X_imp = tmp_X,
                                  Z_imp = tmp_Z,
                                  clID = data_before[, fe$clID_varname],
                                  heap = 0,
                                  nitt = nitt,
                                  thin = thin,
                                  burnin = burnin)

      }

    }


    if(tmp_type == "roundedcont"){

      #yet, we don't have a multilevel imputation for rounded incomes
      imp <- imp_roundedcont(y_imp = data_before[, l2],
                             X_imp = tmp_X)

    }

    if(tmp_type == "interval"){

      #yet, we don't have a multilevel imputation for interval data
      imp <- imp_interval(y_imp = data_before[, l2],
                          X_imp = tmp_X)

    }


    if(tmp_type == "count"){

      if(use_single_level){

        #print("Currently the single level model.")
        imp <- imp_count_single(y_imp = data_before[, l2],
                                X_imp = tmp_X)


      }else{

        #print("Currently the multilevel model.")
        imp <- imp_count_multi(y_imp = data_before[, l2],
                               X_imp = tmp_X,
                               Z_imp = tmp_Z,
                               clID = data_before[, fe$clID_varname],
                               nitt = nitt,
                               thin = thin,
                               burnin =  burnin)

      }

    }

    if(tmp_type == "categorical"){

      # if the cluster ID is to be imputed, obviously we cannot run a multilevel model.
      if(use_single_level || l2 == fe$clID_varname){

        #print("Currently the single level model.")
        imp <- imp_cat_single(y_imp = data_before[, l2],
                              X_imp = tmp_X)

      }else{

        #print("Currently the multilevel model.")
        imp <- imp_cat_multi(y_imp = as.factor(data_before[, l2]),
                             X_imp = tmp_X,
                             Z_imp = tmp_Z,
                             clID = data_before[, fe$clID_varname],
                             nitt = nitt,
                             thin = thin,
                             burnin = burnin)
      }
    }

    if(tmp_type == "ordered_categorical"){

      if(use_single_level){

        #print("Currently the single level model.")
        imp <- imp_orderedcat_single(y_imp = data_before[, l2],
                                     X_imp = tmp_X)

      }else{

        #print("Currently the multilevel model.")
        imp <- imp_orderedcat_multi(y_imp = as.factor(data_before[, l2]),
                                    X_imp = tmp_X,
                                    Z_imp = tmp_Z,
                                    clID = data_before[, fe$clID_varname],
                                    nitt = nitt,
                                    thin = thin,
                                    burnin = burnin)
      }
    }

    if(tmp_type == "intercept"){
      imp <- sample_imp(data_before[, l2])
    }
    data_before[, l2] <- imp
  }
  data_after <- data_before[, original_ordering, drop = FALSE]
  return(data_after)
}
