

#' Sample imputation.
#'
#' Function to sample values in a variable from other (observed) values in this variable.
#' So this imputation doesn't use further covariates.
#' @param variable A vector of size \code{n} with missing values.
#' @return A vector of size \code{n} without missing values.
#' @examples
#' set.seed(123)
#' sample_imp(c(1, NA, 3, NA, 5))
#' @export
sample_imp <- function(variable){

  if(all(is.na(variable))) stop("Variable consists only of NAs.")
  ret <- variable
  ret[is.na(ret)] <- sample(size = sum(is.na(variable)),
                             variable[!is.na(variable)], replace = TRUE)
  return(ret)
}


#' Get the type of variables.
#'
#' Function checks wether a variable is: ...
#' \itemize{
#'  \item continuous (numeric values),
#'  \item semicontinuous (numeric values with more than 5\% of them are 0),
#'  \item rounded continuous (if continuous values are rounded to the closest multiple of 5, 10, 50, 100, 500 or 1000.
#'  We see this to be the case if more than 50\% of the observations are devisable by 5)
#'  \item count data (if all values are integers).
#'  \item an intercept (the same value for all observations),
#'  \item binary (two different values - like 0s and 1s or "m" and "f"),
#'  \item categorical (the variable is a factor or has more than 3 different values)
#'  \item orderd categorical (the categorical variable is ordered.)
#'}
#'
#' @param variable A variable (vector) from your data set.
#' @return A character denoting the type of \code{variable}.
#' @examples get_type(iris$Sepal.Length); get_type(iris$Species)
#' @export
get_type <- function(variable){


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
    type <- "interval"
    #if the precise part of an interval variable is rounded,
    #the whole variable is considered to be a rounded continous variable

    if(sum(decompose_interval(variable)$precise %% 5 == 0, na.rm = TRUE)/
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
      #I count this variable as semi-continious.

      if(sum(variable == 0, na.rm = TRUE)/
         length(variable[!is.na(variable)]) > 0.05){
        type <- "semicont"
      }

      # if more than 50% of the data are divisable by 5, they are considered to be rounded
      # continous
      if(sum(variable %% 5 == 0, na.rm = TRUE)/
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
#' For example, if a continouos variable as only two observations left, then get_type
#' interprete this as a binary variable and not a continouos.
#' @param data the data.frame also passed to \code{hmi}.
#' @return a list with suggested types. Each list element has the name of a variable
#' in the data.frame. The elements contain a single character denoting the type of the variable.
#' See \code{get_type} for details about the variable types.
#' @export
list_of_types_maker <- function(data){

  if(!is.data.frame(data)) stop("We need the data in the data.frame format!")
  ret <- lapply(data, get_type)

}


#' Averages the results of the imputation function \code{hmi}.
#'
#' This function applies the analysis the user is interested in, on all different imputed dataset.
#' Then the results are pooled by simply averaging the results. So the user has to make sure that
#' his analysis produces results with a meaningful average. And furthermore has to accept that no
#' variance is calculated for these parameters.
#' @param mids A \code{mids} (multiple imputed data set) object.
#' Either from the \code{hmi} imputation function or \code{mice}.
#' @param analysis_function A user generated function that contains the model and the model parameters
#' he is interested in. See examples.
#' @return A vector with all averaged results.
#' @examples
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
#' hmi_imp <- hmi(data = test, model_formula = my.formula)
#' hmi_pool(mids = hmi_imp, analysis_function = my_analysis)
#' #if you are interested in fixed effects only, consider using \code{pool} from \code{mice}.
#' pool(with(data = hmi_imp, expr = lmer(Reaction ~ Days + (1 + Days | Subject))))
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
#' (regardles whether they are observed precisly or presumably rounded).
#' @param para This is the vector \eqn{Psi} of parameters
#' (see p. 62 in Drechsler, Kiesl & Speidel, 2015).
#' With respect to them, the value returned by negloglik2 shall be
#' maximized.\cr
#' The starting values are c(kstart, betastart2, gammastart, sigmastart2)
#' (the 6 treshholds (or "cutting points") for the latent variable behind the rounding degree,
#' the regression parameters explaining the logged income,
#' the regression parameters explaining the rounding degree
#' and the variance parameter).
#' @param X the data.frame of covariates.
#' @param y_in_negloglik the target variable (a vector).
#' @param lower the lower bound of an interval variable.
#' @param upper the upper bound of an interval variable.
#' @param my_p This vector is the indicator of the (highes possible) rounding degree for an observation.
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

  # the first 6 parameters are the tresholds for the rounding degree
  # below k0: no rounding, between k0 and k1: rounding degree 5 etc.
  # The other degrees are 10, 50, 100, 500 and 1000
  k0 <- para[1]
  k1 <- para[2]
  k2 <- para[3]
  k3 <- para[4]
  k4 <- para[5]
  k5 <- para[6]

  # the regression coefficients beta defining mean2 - the expected value of log(Y)
  # (see eq (2) in Drechsler, Kiesl & Speidel, 2015).
  # They also appear in mean1 (the expected value of G).
  reg_coefs <- matrix(para[7:(length(para) - 2)], ncol = 1)

  gamma1 <- para[length(para) - 1] #the gamma1 from eq (3)

  sigmax <- para[length(para)] # the variance for log(Y) (see eq (3) in Drechsler, Kiesl & Speidel, 2015)

  if (k1 < k0) return(1e50)        # make sure that threshholds
  if (k2 < k1) return(1e50)        # are in increasing order
  if (k3 < k2) return(1e50)        #
  if (k4 < k3) return(1e50)        #
  if (k5 < k4) return(1e50)        #

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
  intervall_obs <- !is.na(lower)
  sum_likelihood_imprecise_obs <- sum(log(contributions4intervalls(lower = lower[intervall_obs],
                                                                   upper = upper[intervall_obs],
                                                                   mymean = mean2[intervall_obs],
                                                                   mysd = sigma2)))

  # a1 are a7 the "components" of the log-liklihood contributions
  #
  # If the income is not divisable by 5, only a1 is relevant
  # if the income is divisable by 5, but not by 10, a1 + a2 are relevant
  # etc.

  #get the value of the standard bivariate normal cumulative density funtion
  #the better k0, k1, ..., sigma and rho are found, the better a1, a2, ...

  a1 <- pbivnormX(x =  c(k0 - mean1[!missind_negloglik])/sigma1,
                  y =  ((log(y_obs + 0.5) - mean_y_precise)/sd_y_precise - mean2[!missind_negloglik])/sigma2,
                  rho = rho
  ) -
    pbivnormX(x =  c(k0 - mean1[!missind_negloglik])/sigma1,
              y =  ((log(pmax(y_obs - 0.5, 0)) - mean_y_precise)/sd_y_precise - mean2[!missind_negloglik])/sigma2,
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
                   half_divisor = 2.5)

  a3 <- components(ki = k1, kj = k2,
                   mean1_obs = mean1[!missind_negloglik],
                   mean2_obs = mean2[!missind_negloglik],
                   sigma1 = sigma1,
                   sigma2 = sigma2,
                   rho = rho,
                   y_obs = y_obs,
                   mean_y_precise = mean_y_precise,
                   sd_y_precise = sd_y_precise,
                   half_divisor = 5)

  a4 <- components(ki = k2, kj = k3,
                   mean1_obs = mean1[!missind_negloglik],
                   mean2_obs = mean2[!missind_negloglik],
                   sigma1 = sigma1,
                   sigma2 = sigma2,
                   rho = rho,
                   y_obs = y_obs,
                   mean_y_precise = mean_y_precise,
                   sd_y_precise = sd_y_precise,
                   half_divisor = 25)

  a5 <- components(ki = k3, kj = k4,
                   mean1_obs = mean1[!missind_negloglik],
                   mean2_obs = mean2[!missind_negloglik],
                   sigma1 = sigma1,
                   sigma2 = sigma2,
                   rho = rho,
                   y_obs = y_obs,
                   mean_y_precise= mean_y_precise,
                   sd_y_precise= sd_y_precise,
                   half_divisor = 50)

  a6 <- components(ki = k4, kj = k5,
                   mean1_obs = mean1[!missind_negloglik],
                   mean2_obs = mean2[!missind_negloglik],
                   sigma1 = sigma1,
                   sigma2 = sigma2,
                   rho = rho,
                   y_obs = y_obs,
                   mean_y_precise = mean_y_precise,
                   sd_y_precise = sd_y_precise,
                   half_divisor = 250)


  a7 <- stats::pnorm((log(y_obs + 500) - mean_y_precise)/sd_y_precise, mean2[!missind_negloglik], sigma2) -

    pbivnormX( x =  c(k5 - mean1[!missind_negloglik])/sigma1,
               y =  ((log(y_obs + 500) - mean_y_precise)/sd_y_precise - mean2[!missind_negloglik])/sigma2,
               rho = rho
    ) -

    stats::pnorm((log(pmax(y_obs - 500, 0)) - mean_y_precise)/sd_y_precise, mean2[!missind_negloglik], sigma2)  +

    pbivnormX( x =  c(k5 - mean1[!missind_negloglik])/sigma1,
               y =  ((log(pmax(y_obs - 500, 0)) - mean_y_precise)/sd_y_precise - mean2[!missind_negloglik])/sigma2,
               rho = rho
    )

  a2 <- pmax(0, a2) # replacing negative values by 0.
  a3 <- pmax(0, a3)
  a4 <- pmax(0, a4)
  a5 <- pmax(0, a5)
  a6 <- pmax(0, a6)
  a7 <- pmax(0, a7)

  result_try <- cbind(a1 , a2 * (my_p[!missind_negloglik] >= 1) , a3 * (my_p[!missind_negloglik] >= 2),
                      a4 * (my_p[!missind_negloglik] >= 3) , a5 * (my_p[!missind_negloglik] >= 4),
                      a6 * (my_p[!missind_negloglik] >= 5), a7 * (my_p[!missind_negloglik] == 6))

  # Sum up the likelihood contributions (see eq (5) in Drechsler, Kiesl & Speidel, 2015)
  result <- log(rowSums(result_try, na.rm = TRUE))

  ret <- sum(result) + sum_likelihood_imprecise_obs
  return(-ret) # notice the minus
}

#' calculate the likelihood contribution of the intervall data only
#'
#' calculate the likelihood contribution of the intervall data only
#' @param para This is the vector \eqn{Psi} of parameters
#' (see p. 62 in Drechsler, Kiesl & Speidel, 2015).
#' With respect to them, the value returned by negloglik2 shall be
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
  intervall_obs <- !is.na(lower)
  sum_likelihood_imprecise_obs <- sum(log(contributions4intervalls(lower = lower[intervall_obs],
                                                                   upper = upper[intervall_obs],
                                                                   mymean = mean2[intervall_obs],
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
#' @param mysd a numeric for the standard devitation of the normal distribution
#' (which is one of the parameters trying to be optimized
#' so that the likelihood becomes maximized)
#' @return a vector giving the likelihood contributions of the interval data.
contributions4intervalls <- function(lower, upper, mymean, mysd){

  ret <- stats::pnorm(upper, mean = mymean, sd = mysd) -
    stats::pnorm(lower, mean = mymean, sd = mysd)

  return(ret)
}

#' calculate probabilities from the cumulative distribution function of a standard bivariate normal distribution
#'
#' A modfied version of pbivnorm() from package \code{pbivnorm}.
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


#' Function to get the likelihood contribution of different roudinding degrees
#'
#' It is needed in the imputation routine for rounded income.
#' See equation (5) in Drechsler, Kiesel & Speidel (2015)
#' @param ki An integer for the i-th treshold (or "cutpoint")
#' @param kj An integer for the j-th treshold (or "cutpoint") (ki < kj)
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

  upper_inner <- ((log(y_obs + half_divisor) - mean_y_precise)/sd_y_precise - mean2_obs)/sigma2
  lower_inner <- ((log(pmax(y_obs - half_divisor, 0)) - mean_y_precise)/sd_y_precise - mean2_obs)/sigma2
  lower_outer <- c(ki - mean1_obs)/sigma1
  upper_outer <- c(kj - mean1_obs)/sigma1

  ret <- doubleintegral(lower_inner = lower_inner, upper_inner = upper_inner,
                        lower_outer = lower_outer, upper_outer = upper_outer,
                        cdf = pbivnormX, rho = rho)

  return(ret)
}

#' Function to calculate double integrals
#'
#' This function is primarily build to make the funtion \code{components} more neat.
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


#' Function to extract the different elements of a formula
#'
#' The function searches for the target variable, fixed effects variables,
#' if there is a cluster ID: this and the random effects variables.\cr
#' The names of the fixed and random intercepts variable (if existent) are explicitly labelled.
#' @param model_formula A formula (from class \code{formula})
#' @param constant_variables A Boolean-vector of length equal to the number of columns in the data set
#' specifying whether a variable is a constant variable (eg. an intercept variable) or not.
#' @param variable_names_in_data A character-vector with the column names of the data set.
#' @return A list with the names of the target variable, the intercept variable,
#' the fixed and random effects covariates and the cluster id variable.\cr
#' If some of them doesn't exist, they get the value "".
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
    # get a good name for this intercept variable and make shure that an
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

    intercept_varname <- ""
    fixedeffects_varname <- variable_names_in_data
    randomeffects_varname <- ""
    clID_varname <- ""
  }

  # Note: there can be only one intercept variable in the data set.

  ret <- list(target_varname = target_varname,
              intercept_varname = intercept_varname,
              fixedeffects_varname = fixedeffects_varname,
              randomeffects_varname = randomeffects_varname,
              clID_varname = clID_varname)
  return(ret)
}

#' Function to transform numeric (or character) vectors into an interval object
#'
#' Function to transform numeric (or character) vectors into an interval object
#' @param lower a vector with the lower values of a variable
#' @param upper a vector with the upper values of a variable
#' @return a character vector with the lower and upper bound values.
#' @export
as_interval <- function(lower, upper){
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


#' Function to check wheter an object is an interval
#'
#' If there are numerics seperated by semicolons (;), this is considered to be an interval.
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
  reg_low <- "^((-){0,1}[[:digit:]]{1,}(\\.[[:digit:]]{1,}){0,1}|-Inf)"

  #the end of the interval is one or more digits optionally preceded by an "-"
  # optionally followed by one period (or "full stop") and one or more digits
  #or the interval ends with "Inf"
  reg_up <- "((-){0,1}[[:digit:]]{1,}(\\.[[:digit:]]{1,}){0,1}|Inf)$"
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
as_data_frame_interval <- function(x, ...){
  tmp <- as.data.frame(as.character(x), ...)
  colnames(tmp) <- "interval"
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
#' @param interval an interval object of length n
#' (if it is something else, it is returned unchanged)
#' @return a data.frame with the precise observations, the lower and the upper bounds
#' as covariates.
decompose_interval <- function(interval){
  if(!is_interval(interval)){
    return(data.frame(precise = interval, lower = -Inf, upper = Inf))
  }
  tmp <- split_interval(interval)
  dists <- tmp[, 2] - tmp[, 1]
  precise <- ifelse(abs(dists) < 1e-20, tmp[, 1], NA)
  imprecise <- ifelse(abs(dists) > 1e-20, interval, NA)
  tmp <- split_interval(imprecise)
  return(data.frame(precise, lower = tmp[, 1], upper = tmp[, 2]))
}

#' Get standard NAs
#'
#' This function replaces observations with "-Inf;Inf" with the standard NAs (therefore 'sna')
#' @param x an object. Designed for \code{interval}-objects.
#' (if it is something else, it is returned unchanged)
#' @return a n times 2 matrix. The first column is the lower bound, the second the upper bound.
sna_interval <- function(x){

  #if it is no interval object, just return the object
  if(!is_interval(x)) return(x)

  #NA values have to be set temporarily to "NA;NA"
  x[x == "-Inf;Inf"] <- NA
  return(x)
}

#' Transform interval variables to a interval data frame
#'
#' This function is the path from this hmi package to the linLIR package (Wiencierz, 2012).
#' @param interval an \code{interval}
#' @return an interval data frame (idf-object) with one variable (having a lower and a upper bound).
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
    ret[, i] <- as_interval(lower = idf[[i]][, 1], upper = idf[[i]][, 2])
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
'+.interval' <- function(interval, x){
  tmp <- split_interval(interval) + split_interval(x)
  as_interval(lower = tmp[, 1], upper = tmp[, 2])
}


#' Subraction function
#'
#' Function to subtract single elements or vectors (of correct dimension) from the interval object
#' @param interval an object from class interval
#' @param x an single element or vector to subtract to the interval object
#' @return an interval object
#' @export
'-.interval' <- function(interval, x){
  tmp <- split_interval(interval) - split_interval(x)
  as_interval(lower = tmp[, 1], upper = tmp[, 2])
}

#' Multiplication function
#'
#' Function to multiply single elements or vectors (of correct dimension) to the interval object
#' @param interval an object from class interval
#' @param x an single element or vector for multplication
#' @return an interval object
#' @export
'*.interval' <- function(interval, x){
  tmp <- split_interval(interval) * split_interval(x)
  as_interval(lower = tmp[, 1], upper = tmp[, 2])
}

#' Dividing function
#'
#' Function to divide single elements or vectors (of correct dimension) to the interval object
#' @param interval an object from class interval
#' @param x an single element or vector for division
#' @return an interval object
#' @export
'/.interval' <- function(interval, x){
  tmp <- split_interval(interval) / split_interval(x)
  as_interval(lower = tmp[, 1], upper = tmp[, 2])
}


