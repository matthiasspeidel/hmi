#' Sample imputation.
#'
#' Function to sample values in a variable from other (observed) values in this variable.
#' So this imputation does not use further covariates.
#' @param variable A vector of size \code{n} with missing values.
#' @return A list with a n times 1 data.frame without missing values and
#'  a list with the chains of the Gibbs-samples for the fixed effects and variance parameters.
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
#' @param rounding_degrees A numeric vector with the presumed rounding degrees.
#' @return A character denoting the type of \code{variable}.
#' @examples get_type(iris$Sepal.Length); get_type(iris$Species)
#' @export
get_type <- function(variable,
                     rounding_degrees = NULL){

  if(is.matrix(variable) && ncol(variable) == 1){
    variable <- variable[, 1]
  }
  type <- "unknown"
  if(all(is.na(variable))) return(type)
  if(!(is.vector(rounding_degrees) | is.null(rounding_degrees))){
    stop("rounding_degrees passed to get_type() have to be a vector or NULL.")
  }



  if(length(unique(variable[!is.na(variable)])) == 1){
    type <- "intercept"
    return(type)
  }

  if(length(unique(variable[!is.na(variable)])) == 2){
    type <- "binary"

    return(type)
  }


  if(is.integer(variable)){
    type <- "count"

    return(type)
  }


  if(is.vector(variable) || is.factor(variable) || is.array(variable) || is.numeric(variable)){

    #checking the length of the table to be larger than two is not sufficient
    #to get a real categoriacal variable, because for a (real) continious variable
    #the length of the table is equal to the length of the variable
    #(what will be in general larger than two...)
    #Therefore it will be checked if the table has more than two entries AND
    #wheter the varriable is a factor.
    if(is.factor(variable)){
      type <- "categorical"
      if(is.ordered(variable)){
        type <- "ordered_categorical"
      }

      return(type)
    }

    #If no default rounding_degrees were given, suggest_rounding_degrees suggests them
    if(is.null(rounding_degrees)){
      rounding_degrees <- suggest_rounding_degrees(variable)
    }

    #... if they are still NULL, then c(1, 10, 100, 1000) is used as a default
    if(is.null(rounding_degrees)){
      rounding_degrees <- c(1, 10, 100, 1000)
    }

    if(is.numeric(variable)){
      type <- "cont"

      # if there is no difference between a variable and their as.integer()-result,
      # it is considered to be an integer...
      if(max(abs(variable - as.integer(variable)), na.rm = TRUE) == 0){
        type <- "count"

        #... but not if too many categories are available...
        if(length(table(variable)) > 20){
          type <- "cont"
        }

        #... or it is a rounded continuous variable:
        if(sum(apply(outer(decompose_interval(interval = variable)[, "precise"], rounding_degrees[-1], '%%') == 0,
                     1, any), na.rm = TRUE)/
           length(variable[!is.na(variable)]) > 0.5){
          type <- "roundedcont"
        }

        return(type)
      }

      #if the variable is numeric and more than 10% of the variable values share the same value,
      #We consider this variable as semi-continious.

      if(max(table(variable))/length(variable[!is.na(variable)]) > 0.1){
        type <- "semicont"
      }

      # if more than 50% of the data are divisible by on of the given rounding degrees,
      # they are considered to be rounded continuous
      # Observed 0s shall not be counted as rounded,
      #because otherwise semi-continuous variables might be considered to be rounded-continuous.

      if((sum(apply(outer(decompose_interval(interval = variable)[, "precise"], rounding_degrees[-1], '%%') == 0,
                    1, any), na.rm = TRUE) - sum(variable == 0, na.rm = TRUE))/
        length(variable[!is.na(variable)]) > 0.5){
        type <- "roundedcont"
      }

      if(type == "unknown"){
        type <- "categorical"
      }
      return(type)

    }



  }else{

    if(is_interval(variable)){

      #If no default rounding_degrees were given, suggest_rounding_degrees suggests them
      if(is.null(rounding_degrees)){
        rounding_degrees <- suggest_rounding_degrees(variable)
      }

      #... if they are still NULL, then c(1, 10, 100, 1000) is used as a default
      if(is.null(rounding_degrees)){
        rounding_degrees <- c(1, 10, 100, 1000)
      }

      type <- "interval"
      #if more than 50% of the precise part of an interval variable is rounded,
      #the whole variable is considered to be a rounded continous variable
      tmp <- decompose_interval(interval = variable)[, "precise"]
      if(sum(!is.na(tmp)) == 0) return(type)
      if(sum(apply(outer(tmp, rounding_degrees[-1], '%%') == 0,
                   1, any), na.rm = TRUE)/
         length(tmp[!is.na(tmp)]) > 0.5){
        type <- "roundedcont"
      }
      return(type)
    }

    if(ncol(variable) == 1){
      ret <- get_type(variable[, 1], rounding_degrees = rounding_degrees)
      return(ret)
    }
    if(ncol(variable) > 1){
      if(is.matrix(variable)){
        variable <- as.data.frame(variable)
      }
      ret <- unlist(lapply(variable, get_type, rounding_degrees = rounding_degrees))

      return(ret)
    }
  }
  return(type)
}

#' Get the mode
#'
#' This function calculates the mode (most frequent observation) of a vector.
#' @param x A vector
#' @references Adopted from stackoverflow.com/questions/2547402: "is there a built in function for finding the mode"
#' from user "Ken Williams".
#' @return The mode of x as a numeric value.
Mode <- function(x){
  ux <- unique(x[!is.na(x)])
  return(ux[which.max(tabulate(match(x, ux)))])
}

#' Helps the user to make a list of heaps.
#'
#' In \code{hmi} the user can add a list of heaps. This function gives him a framework
#' with suggestions. Of course the user can make changes by herself/himself afterwards.
#' For example, the function might wrongly classify a variable to be heaped.
#' @param data the data.frame also passed to \code{hmi}.
#' @return a list with suggested heaps. Each list element has the name of a heaped variable
#' in the data.frame. The elements contain a single numeric denoting the heap found for that variable.
#' @export
list_of_spikes_maker <- function(data){
  if(!is.data.frame(data)) stop("We need the data in the data.frame format!")
  ret <- list()
  for(l2 in colnames(data)){
    if(get_type(data[, l2]) == "semicont"){
      ret[[l2]] <- Mode(data[, l2])
    }

  }
  return(ret)
}

#' Helps the user to make a list of types.
#'
#' In \code{hmi} the user can add a list of types. This function gives him a framework
#' with suggestions. Of course the user can make changes by herself/himself afterwards.
#' For example, if a continuous variable as only two observations left, then get_type
#' interpret this as a binary variable and not a continuous.
#' @param data the data.frame also passed to \code{hmi}.
#' @param rounding_degrees A numeric vector with the presumed rounding degrees.
#' @return a list with suggested types. Each list element has the name of a variable
#' in the data.frame. The elements contain a single character denoting the type of the variable.
#' See \code{get_type} for details about the variable types.
#' @export
list_of_types_maker <- function(data, rounding_degrees = NULL){

  if(!is.data.frame(data)) stop("We need the data in the data.frame format!")

  rounding_degrees_tmp <- rounding_degrees

  ret <- list()

  for(l2 in colnames(data)){
    if(is.list(rounding_degrees)){
      rounding_degrees_tmp <- rounding_degrees[[l2]]
    }

    ret[[l2]] <- get_type(data[, l2], rounding_degrees = rounding_degrees_tmp)
  }

  return(ret)
}

#' Function to get all factors
#'
#' Function to get all factors (not limited to prime factors) of an integer.
#' @param x An integer.
#' @return A numeric vector with the factors
#' @references based on stackoverflow.com/questions/6424856 "R Function for returning ALL factors"
#'  answer by Chase
factors <- function(x){
  if(is.na(x)) return(NA)
  if(x %% 1 != 0) return(NA)
  x <- as.integer(x)
  div <- seq_len(abs(x))
  return(div[x %% div == 0L])
}


#' suggesting rounding degrees
#'
#' A function that suggests some rounding degrees of a continuous variable
#' (classically formated or as interval object)
#' @param x A vector or \code{interval} object.
suggest_rounding_degrees <- function(x){

  if(!(is.vector(x) | is_interval(x))){
    return(NULL)
  }

  variable <- decompose_interval(interval = x)[, "precise"]
  tab <- table(unlist(sapply(variable, factors)))

  if(length(tab) == 0){
    return(NULL)
  }else{

    values <- as.numeric(names(tab))

    #for later, only those values shall be considered if the observed number of indiduvials,
    #rounded to this degree, exceeds the expected number at least by two.
    #Example: from 10000 individuals, 2000 are expected to round to 5.
    #If 4000 or more indivudials are observed 4 is considered to be a possible rounding degree
    candidates_values <- values[values * tab/ length(x) > 2]
    if(length(candidates_values) == 0) return(NULL)
    sorted_candidates_values <- sort(candidates_values, decreasing = TRUE)
    #now check the relative frequencies of these rounding degrees...
    tmp <- outer(variable, sorted_candidates_values, '%%') == 0
    colnames(tmp) <- sorted_candidates_values

    rel_rounding_freqs <- colSums(tmp, na.rm = TRUE)/sum(!is.na(variable))
    current_percent <- 0
    for(i in 1:ncol(tmp)){
      #...if a rounding degree does provide less then 10 percent new people rounded to this degree,
      #it is considered to be a neglictible rounding degree
      if(rel_rounding_freqs[i] - current_percent < 0.1){
        rel_rounding_freqs[i] <- NA
      }else{
        #Here we implicity assume, that smaller rounding degrees are a factor of the previous ones
        #(eg. 10 as a factor of 100 and 1000).
        #if multiple bases (eg. 10 and 7) and there follow ups (like 100, 1000, 14, 21)
        #should be found, these bases would beed an own current_percent value.
        #Currently this seems to be an overkill.
        current_percent <- rel_rounding_freqs[i]
      }
    }

    suggested_rounding_degrees <- c(1, sort(sorted_candidates_values[!is.na(rel_rounding_freqs)]))
    if(get_type(variable, rounding_degrees = suggested_rounding_degrees) == "roundedcont"){
      return(suggested_rounding_degrees)
    }else{
      return(NULL)
    }
  }
}


#' Helps the user to make a list of rounding degrees
#'
#' In \code{hmi} the user can add a list of rounding degrees. This function gives him a framework
#' with suggestions. Of course the user can make changes by herself/himself afterwards.
#' For example, the function might wrongly classify a variable to be heaped or selects
#' unwanted rounding degrees.
#' @param data the data.frame also passed to \code{hmi}.
#' @return a list with suggested rounding degrees. Each list element has the name
#' of a rounded continuous variable in the data.frame. The elements contain a
#' numeric vector with the rounding degrees found for that variable.
#' @export
list_of_rounding_degrees_maker <- function(data){
  if(!is.data.frame(data)) stop("We need the data in the data.frame format!")
  ret <- list()

  for(l2 in colnames(data)){

    variable <- decompose_interval(interval = data[, l2])[, "precise"]
    suggested_rounding_degrees <- suggest_rounding_degrees(variable)
    if(get_type(variable, rounding_degrees = suggested_rounding_degrees) == "roundedcont"){
      ret[[l2]] <- suggested_rounding_degrees
    }
  }
  return(ret)
}


#' Helps the user to make a list of rounding formulas for the rounding degress
#'
#' In \code{hmi} the user can add a list of heaps. This function gives him a framework
#' with suggestions. Of course the user can make changes by herself/himself afterwards.
#' For example, the function might wrongly classify a variable to be heaped.
#' @param data the data.frame also passed to \code{hmi}.
#' @param default A default formula used for every rounded variable.
#' @return a list with suggested rounding degree formulas. Each list element has the name
#' of a rounded continuous variable in the data.frame. The elements contain a
#' very general rounding degree formula.
#' @export
list_of_rounding_formulas_maker <- function(data, default = ~ .){
  if(!is.data.frame(data)) stop("We need the data in the data.frame format!")
  ret <- list()

  for(l2 in colnames(data)){
    variable <- decompose_interval(interval = data[, l2])[, "precise"]
    suggested_rounding_degrees <- suggest_rounding_degrees(variable)
    if(get_type(variable, rounding_degrees = suggested_rounding_degrees) == "roundedcont"){
      ret[[l2]] <- stats::formula(default)
    }
  }
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
#' With respect to them, the value returned by negloglik shall be
#' maximized.\cr
#' The starting values are c(kstart, betastart, gamma1start, sigmastart)
#' (the thresholds (or "cutting points") for the latent variable behind the rounding degree,
#' the regression parameters explaining the logged income,
#' the regression parameters explaining the rounding degree
#' and the variance parameter).
#' @param X_in_negloglik The data.frame of covariates explaining Y, the observed target variable.
#' @param PSI_in_negloglik The data.frame of covariates explaining G, the latent rounding tendency.
#' Without the target variable.
#' @param vars_in_psi A vector with the names of the variables that should be used in PSI
#' (the variables explaning the latent rounding tendency G), without the intercept and Y.
#' @param y_precise_stand A vector of the precise (and standardized) observations from the target variable.
#' @param lower_bounds The lower bounds of an interval variable.
#' @param upper_bounds The upper bounds of an interval variable.
#' @param my_g This vector is the indicator of the (highest possible) rounding degree for an observation.
#' This parameter comes directly from the data.
#' @param sd_of_y_precise The scalar with the value equal to the standard deviation of the target variable.
#' @param indicator_precise A boolean Vector indicating whether the value in the original target
#' variable is precise (e.g. 5123 or 5123.643634) or not.
#' @param indicator_imprecise A boolean Vector indicating whether the value in the original target
#' variable is imprecise (e.g. "5120;5130) or not.
#' @param indicator_outliers A boolean Vector indicating whether the value in the precise
#' observations of the original target are outliers (smaller than 0.5\% or
#' larger than 99.5\% of the other precise observations).
#' @param rounding_degrees A numeric vector with the presumed rounding degrees.
#' @references Joerg Drechsler, Hans Kiesl, Matthias Speidel (2015):
#' "MI Double Feature: Multiple Imputation to Address Nonresponse and Rounding Errors in Income Questions",
#' Austrian Journal of Statistics, Vol. 44, No. 2, \url{http://dx.doi.org/10.17713/ajs.v44i2.77}
#' @return An integer equal to the (sum of the) negative log-likelihood contributions (of the observations)
negloglik <- function(para,
                      X_in_negloglik,
                      PSI_in_negloglik,
                      vars_in_psi,
                      y_precise_stand,
                      lower_bounds = NA, upper_bounds = NA,
                      my_g, sd_of_y_precise,
                      indicator_precise,
                      indicator_imprecise,
                      indicator_outliers,
                      rounding_degrees = c(1, 10, 100, 1000)){

  lower <- ifelse(is.na(lower_bounds), -Inf, lower_bounds)
  upper <- ifelse(is.na(upper_bounds), Inf, upper_bounds)

  # the first parameters are the thresholds for the rounding degree
  thresholds <- para[grep("^threshold", names(para))]
  # the regression coefficients beta defining mean2 - the expected value of log(Y)
  # (see eq (2) in Drechsler, Kiesl & Speidel, 2015).
  # They might also appear in mean1 (the expected value of G).
  #The intercept is not part of the maximazition, as due to standardizations of y and x,
  #its value is exactly 1.
  #If there is only an intercept variable in X, there are no coefficients to be estimated.
  if(ncol(X_in_negloglik) == 1){
    reg_coefs <- matrix(1, ncol = 1)
  }else{
    #!!!DEBUG!!! REMOVE NULL (OR REPLACE IT WITH 1, IF WE WANT TO SET THE INTERCEPT TO 1)
    reg_coefs <- matrix(c(NULL, para[grep("^coef_y_on_x", names(para))]), ncol = 1)
  }

  G_coefs <- matrix(para[grep("^coef_g_on_psi", names(para))], ncol = 1)

  gamma1 <-  para[grep("^gamma1", names(para))] #the gamma1 from eq (3)
  if(length(gamma1) == 0){
    gamma1 <- 0
  }
  sigmax <- para[grep("^sigma", names(para))] # the (residual) variance for log(Y) (see eq (3) in Drechsler, Kiesl & Speidel, 2015)

  if(is.unsorted(thresholds)) return(1e50) # make sure that threshholds are increasing
  if(sigmax <= 0) return(1e50) # the residual variance has to be positive

  # mean and covariance of the joint normal of x and g ,
  # following Heitjan & Rubin

  #the mean for G (see eq (2) in Drechsler, Kiesl & Speidel, 2015)
  individual_means_for_g <- gamma1 *
    (as.matrix(X_in_negloglik[indicator_precise, , drop = FALSE][!indicator_outliers, , drop = FALSE]) %*% reg_coefs) +
    as.matrix(PSI_in_negloglik[, vars_in_psi, drop = FALSE][!indicator_outliers, , drop = FALSE]) %*% G_coefs

  #the mean for log(Y) (see eq (2) in Drechsler, Kiesl & Speidel, 2015)
  individual_means_for_y <- as.matrix(X_in_negloglik) %*% reg_coefs

  #the covariance matrix for G and log(Y)  (see eq (3) in Drechsler, Kiesl & Speidel, 2015)
  sigma <- matrix(c(1 + gamma1^2 * sigmax^2,
                    gamma1 * sigmax^2,
                    gamma1 * sigmax^2,
                    sigmax^2), nrow = 2)

  rho <- max(min(stats::cov2cor(sigma)[1, 2], 1), -1)

  sigma_for_g <- sqrt(sigma[1, 1])
  sigma_for_y <- sqrt(sigma[2, 2])

  #get the likelihood contributions of the imprecise observations
  sum_likelihood_imprecise_obs <- sum(log(contributions4intervals(lower_bounds = lower,
                                                                  upper_bounds = upper,
                                                                  mymean = individual_means_for_y[indicator_imprecise],
                                                                  mysd = sigma_for_y)))

  #check for all possible rounding degrees, it is checked, what the likelihood conrtibution
  #of G and Y is. Each possible rounding degree comes with a range of values for G and
  #a range of values for Y.
  # The range of values for G is set up by the thresholds k, the range of values for Y
  # by the interval: Y_observed +- half the length of the rounding degree.
  # To keep the code more simple, here for every observations all possible rounding degrees are
  # evaluated (e.g. for an observation 3200 the likelihood contribution in the case
  # of rounding to the next multiple of 1000 is evaluated); the irrelevant likelihood-contributions
  # are set to 0 in a second step.
  y_tmp <- y_precise_stand[!indicator_outliers]
  likelihood_contributions_heaped <- array(dim = c(length(y_tmp), length(rounding_degrees)))

  thresholds_extended <- c(-Inf, thresholds, Inf)

  for(j in 1:ncol(likelihood_contributions_heaped)){

    # pbivnormX calculates the probabilities of the standard bivariate normal distribution.
    # Therefore the values going in, have to be standardized

    # integration bounds for y
    upper_inner <- (y_tmp + rounding_degrees[j]/2/sd_of_y_precise - individual_means_for_y[indicator_precise][!indicator_outliers])/sigma_for_y
    lower_inner <- (y_tmp - rounding_degrees[j]/2/sd_of_y_precise - individual_means_for_y[indicator_precise][!indicator_outliers])/sigma_for_y

    # integrations bounds for g. The function c() is used to get a vector.
    lower_outer <- c(thresholds_extended[j] - individual_means_for_g)/sigma_for_g
    upper_outer <- c(thresholds_extended[j + 1] - individual_means_for_g)/sigma_for_g

    likelihood_contributions_heaped[, j] <- pmax(0,
                                                 doubleintegral(lower_inner = lower_inner, upper_inner = upper_inner,
                                                                lower_outer = lower_outer, upper_outer = upper_outer,
                                                                cdf = pbivnormX, rho = rho))
  }

  merged <- cbind(likelihood_contributions_heaped, my_g[!indicator_outliers])

  result_try <- apply(merged, 1, function(x) sum(x[1:x[length(x)]]))
  result_try[is.na(result_try)] <- 0
  # Sum up the likelihood contributions (see eq (5) in Drechsler, Kiesl & Speidel, 2015)
  sum_likelihood_heaped_obs <- sum(log(result_try), na.rm = TRUE)

  ret <- sum_likelihood_heaped_obs + sum_likelihood_imprecise_obs
  return(-ret) # Notice the minus.
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
#' @param lower_bounds the lower bound of an interval variable.
#' @param upper_bounds the upper bound of an interval variable.
#' @references Joerg Drechsler, Hans Kiesl, Matthias Speidel (2015):
#' "MI Double Feature: Multiple Imputation to Address Nonresponse and Rounding Errors in Income Questions",
#' Austrian Journal of Statistics, Vol. 44, No. 2, \url{http://dx.doi.org/10.17713/ajs.v44i2.77}
#' @return An integer equal to the (sum of the) negative log-likelihood contributions (of the observations)
negloglik2_intervalsonly <- function(para, X, lower_bounds, upper_bounds){

  # the regression coefficients beta defining mean2 - the expected value of log(Y)
  # (see eq (2) in Drechsler, Kiesl & Speidel, 2015).
  reg_coefs <- matrix(para[1:(length(para) - 1)], ncol = 1)

  sigmax <- para[length(para)]

  #ascertain that a degenerated variance will not be considered as optimum
  if(sigmax <= 0) return(1e50)
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
  interval_obs <- !is.na(lower_bounds)
  tmp <- contributions4intervals(lower_bounds = lower_bounds[interval_obs],
                                 upper_bounds = upper_bounds[interval_obs],
                                 mymean = mean2[interval_obs],
                                 mysd = sigma2)
  if(any(tmp == 0)) return(1e50)
  sum_likelihood_imprecise_obs <- sum(log(tmp))
  return(-sum_likelihood_imprecise_obs) # notice the minus
}

#' get the likelihood contributions of interval data
#'
#' This function calculates the likelihood contributions of interval data
#' @param lower_bounds a vector with the lower bounds of an interval covariate.
#' @param upper_bounds a vector with the upper bounds of an interval covariate.
#' @param mymean a numeric for the expected value of the normal distribution
#' (which is one of the parameters trying to be optimized
#' so that the likelihood becomes maximized)
#' @param mysd a numeric for the standard deviation of the normal distribution
#' (which is one of the parameters trying to be optimized
#' so that the likelihood becomes maximized)
#' @return a vector giving the likelihood contributions of the interval data.
contributions4intervals <- function(lower_bounds, upper_bounds, mymean, mysd){

  ret <- stats::pnorm(upper_bounds, mean = mymean, sd = mysd) -
    stats::pnorm(lower_bounds, mean = mymean, sd = mysd)

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
                            mean = elements[c(2, 5)],
                            sigma = Sigma,
                            lower = elements[c(1, 4)],
                            upper = elements[c(3, 6)],
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
#' @param data The data.frame the formula belongs to.
#' @return A list with the names of the target variable, the intercept variable,
#' the fixed and random effects covariates (which includes the name of the target variable),
#' the variables with interactions and the cluster id variable.\cr
#' If some of them don't exist, they get the value "".
#' @export
extract_varnames <- function(model_formula = NULL, constant_variables,
                             variable_names_in_data, data){


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

    variable_names_in_formula <- all.vars(stats::terms(model_formula, data = data))

    #looking for possible interactions
    interactions <- grep(":", x = attr(stats::terms(model_formula, data = data), "term.labels"),
                         value = TRUE)
    if(length(interactions) == 0){
      interactions <- ""
    }

    # ----------Target variable
    target_varname_full <- all.vars(stats::update(model_formula, . ~ 1))

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

    fixedeffects_varname <- all.vars(stats::delete.response(stats::terms(lme4::nobars(model_formula), data = data)))
    # --check if intercept is present
    fixed_intercept_exists <- attributes(stats::delete.response(stats::terms(lme4::nobars(model_formula), data = data)))$intercept == 1

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
  if(is.null(target_varname) || target_varname == "."){
    target_varname <- NULL
  }
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
#' The specification of an intercept by calling a 1-column (e.g. "int")
#' is not counted towards the existence of an intercept.
#' Contradictory inputs like "~ 1 + 0 + X1 + ..." or "~ -1 + 1 + X1 + ..."
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
#' The specification of an intercept by calling a 1-column (e.g. "int")
#' is not counted towards the existence of an intercept.
#' Contradictory inputs like "~ 1 + 0 + X1 + ..." or "~ -1 + 1 + X1 + ..."
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

#' Function to transform objects into an interval object
#'
#' Function to transform numeric (or character) vectors or n times 2 matrices into an interval object
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
    if(is.matrix(x) && ncol(x) == 2){
      return(generate_interval(lower = x[, 1], upper = x[, 2]))
    }else{
      x <- paste(x, ";", x, sep = "")
    }

  }
  x[is.na(x)] <- c("NA;NA")
  raw <- unlist(strsplit(x, ";"))
  raw[raw == "NA"] <- NA
  lower <- as.numeric(raw[seq(1, 2*length(x), by = 2)])
  upper <- as.numeric(raw[seq(2, 2*length(x), by = 2)])
  ret <- generate_interval(lower, upper)
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

  if("formula" %in% class(x)) return(FALSE)
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
  tmp <- as.data.frame(rep(NA, length(x)), ...)
  colnames(tmp) <- colnames(x)
  tmp[, 1] <- x
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
#' @param interval A vector, factor or optimally an \code{interval} object of length n
#' (if it is something else, it is returned unchanged)
#' @return A matrix with 5 columns.
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
#' c("2500;2500", "3000;4000", "-Inf;0", NA) will lead to c(2500, 4000, 0, Inf)
decompose_interval <- function(interval){

  if(!(is.vector(interval) | is.factor(interval) | is_interval(interval))){
    stop("interval has to be a vector, factor or interval.")
  }

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
  dists[is.na(dists)] <- Inf
  precise <- ifelse(abs(dists) < 1e-20, tmp[, 1], NA)
  imprecise <- as.interval(ifelse(abs(dists) >= 1e-20 | is.na(dists), interval, NA))
  tmp_imprecise <- split_interval(imprecise)

  ret <- matrix(c(precise, tmp_imprecise[, 1],
                  tmp_imprecise[, 2],
                  ifelse(is.na(tmp[, 1]), -Inf, tmp[, 1]),
                  ifelse(is.na(tmp[, 2]), Inf, tmp[, 2])), ncol = 5)
  colnames(ret) <- c("precise", "lower_imprecise", "upper_imprecise", "lower_general", "upper_general")
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

#' Plotting interval variables
#'
#' Function to plot interval variables by rectangles. The bounds of the rectangles are given by the
#' lower and upper bounds of the interval variables. To avoid precise observations to have a line-width
#' of 0, small values are added to the upper and lower bounds what guarantees the rectangles (or lines or points)
#' to be easily visible in the plot.
#' @param x In its most save way, \code{x} is an object from class \code{interval}
#' and jointly used with a second \code{interval} object \code{y}. If no \code{y} is given,
#' the values of \code{x} are just plotted in order of appearance
#' (cf. \code{plot(iris$Sepal.Length)}).
#' \code{x} can also be a \code{formula} with two variables found in \code{data}.
#' @param y If used jointly with \code{x}, it has to be a numeric vector or an \code{interval} object.
#' @param data If \code{x} is a \code{fomula}, it has to be a data.frame or matrix with column names
#' fitting to the two variables named in the formula.
#' @param col The color of the rectangles.
#' @param xlab A title for the x axis: see \code{title}.
#' @param ylab A title for the y axis: see \code{title}.
#' @param xlim Numeric vectors of length 2, giving the x coordinate ranges.
#' @param ylim Numeric vectors of length 2, giving the y coordinate ranges.
#' @param sort A character specifiying how the values should be sorted if only one variable is to be plotted.
#' By default they are sorted according to their position in the data set.
#' Currently, the only option to chose (\code{sort = "mostprecise_increasing"}) is to sort them by their length of the interval they represent,
#' and within equal lengths increasing with the lower bound.
#' @param ... graphical parameters such as \code{main}.
#' @examples
#' \dontrun{
#' #Works like plot:
#' plot.interval(Sepal.Length ~ Sepal.Width, data = iris)
#' #But designed to plot interval objects:
#' plot.interval(x = artificial$age, y = artificial$income)
#' }
#' @export
plot.interval <- function(x = NULL, y = NULL, data = NULL, col = "black",
                          xlab = NULL, ylab = NULL,
                          xlim = NULL, ylim = NULL,
                          sort = NULL, ...){
  myylab <- "y"
  myxlab <- "x"
  if(is.null(x) && is.null(y)){
    y <- data[, 1]
    x <- data[, 2]
    myylab <- colnames(data)[ , 1]
    myxlab <- colnames(data)[ , 2]
  }

  if((is_interval(x) | is.numeric(x)) && is.null(y)){
    y <- x
    x <- 1:length(y)
    myylab <- "values"
    myxlab <- "Index"
  }

  if(class(x) == "formula"){
    yname <- as.character(x)[2]
    xname <- as.character(x)[3]
    y <- data[, yname]
    x <- data[, xname]
    myylab <- yname
    myxlab <- xname
  }

  if(is_interval(y)){
    y_split <- split_interval(y)
  }else{
    y_split <- cbind(y, y)
  }

  if(is_interval(x)){
    x_split <- split_interval(x)
    x_eps <- (xrange[2] - xrange[1])/500
  }else{
    x_split <- cbind(x, x)
    x_eps <- 0.3
  }

  xrange <- range(x_split, na.rm = TRUE, finite = TRUE)
  yrange <- range(y_split, na.rm = TRUE, finite = TRUE)

  #rect(xleft = -1, ybottom = 1, xright = -0.5, ytop = 1.5, density = 1000)#add ,...
  #make sure points and lines can be seen
  y_eps <- (yrange[2] - yrange[1])/500
  y_split <- y_split + rep(y_eps*c(-1, 1), each = nrow(y_split))


  x_split <- x_split + rep(x_eps*c(-1, 1), each = nrow(x_split))

  #replaces NAs by +- Inf
  x_split[is.na(x_split[, 1]), 1] <- -Inf
  x_split[is.na(x_split[, 2]), 2] <- Inf

  y_split[is.na(y_split[, 1]), 1] <- -Inf
  y_split[is.na(y_split[, 2]), 2] <- Inf

  #replace infinite values with values outside the range
  y_split[is.infinite(y_split[, 1]), 1] <- yrange[1] - (yrange[2] - yrange[1])
  y_split[is.infinite(y_split[, 2]), 2] <- yrange[2] + (yrange[2] - yrange[1])

  x_split[is.infinite(x_split[, 1]), 1] <- xrange[1] - (xrange[2] - xrange[1])
  x_split[is.infinite(x_split[, 2]), 2] <- xrange[2] + (xrange[2] - xrange[1])
  if(!is.null(xlab)) myxlab <- xlab
  if(!is.null(ylab)) myylab <- ylab
  if(!is.null(xlim)) xrange <- xlim
  if(!is.null(ylim)) yrange <- ylim

  if(!is.null(sort)){
    if(sort == "mostprecise_increasing"){
      interval_lengths <- decompose_interval(y)[, "upper_general"] -
        decompose_interval(y)[, "lower_general"]
      #interval_lengths <- y_split[, 2] - y_split[, 1]
      #step1 <- y_split[order(interval_lengths),, drop = FALSE]
      y_split <- y_split[order(interval_lengths, y_split[, 1]),]
      if(is.null(xlab)) myxlab <- "Index of sorted values"
    }
  }



  graphics::plot(0, pch = '', ylab = myylab, xlab = myxlab,
       xlim = xrange,
       ylim = yrange, ...)
  graphics::rect(xleft = x_split[, 1], ybottom = y_split[, 1],
                 xright = x_split[, 2], ytop = y_split[, 2],
       border = NA, col = col, ...)
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
#' @param rounding_degrees A numeric vector with the presumed rounding degrees.
#' @return A n times p data.frame with the standardized versions of the numeric variables.
#' @export
stand <- function(X, rounding_degrees = NULL){
  if(!is.data.frame(X)) stop("X has to be a data.frame.")
  types <- array(NA, dim = ncol(X))
  for(i in 1:length(types)){
    types[i] <- get_type(X[, i], rounding_degrees = rounding_degrees)
  }
  need_stand_X <- types %in% c("cont", "count")
  X_stand <- X
  X_stand[, need_stand_X] <- scale(X[, need_stand_X])
  return(X_stand)
}

#' cleanup data.frames
#'
#' Function to exclude variables that have too many levels
#' (as they may cause numerical problems) and to change binary factors to 0-1 coding
#' (as such factors might generate linear dependent variables).
#' @param X A n times p data.frame with p fixed (or random) effects variables.
#' @param k An integer defining the allowed maximum of levels in a factor covariate.
#' @return A n times (p-r) data.frame, with r being the number of variables with too many factors.
#' @export
cleanup <- function(X, k = 10){

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

  #Make any binary data coded in 0 and 1.
  binary <- types == "binary"
  for(l in which(binary)){
    #make every binaray variable a factor, and then numeric, which should results in
    #a 1 and 2 coding.
    tmp <- as.numeric(factor(X[, l]))
    #With the following formula variables from 1 and 2 coding will be transformed
    #into 0-1 coding.
    #The formula is kept general for the unlikely event that the line above didn't
    #gave a 1-2 but r-s coding.
    #If r < s then the formula recodes r to (r-r)/(s-r) = 0 and s to (s-r)/(s-r) = 1.
    #If r > s then r will be 1 and s = 0.
    X[, l] <- (tmp - min(tmp))/(max(tmp) - min(tmp))
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
#' @param list_of_types a list where each list element has the name of a variable
#' in the data.frame. The elements have to contain a single character denoting the type of the variable.
#' See \code{get_type} for details about the variable types.
#' With the function \code{list_of_types_maker}, the user can get the framework for this object.
#' In most scenarios this is should not be necessary.
#' One example where it might be necessary is when only two observations
#' of a continuous variable are left - because in this case \code{get_type}
#' interpret is variable to be binary. Wrong is it in no case.
#' @param NA_locator A n x p matrix localizing the missing values in the original
#' dataset. The elements are TRUE if the original data are missing and FALSE if the
#' are observed.
#' @param nitt An integer defining number of MCMC iterations (see \code{MCMCglmm}).
#' @param burnin burnin A numeric value between 0 and 1 for the desired percentage of
#' Gibbs samples that shall be regarded as burnin.
#' @param thin An integer to set the thinning interval range. If thin = 1,
#' every iteration of the Gibbs-sampling chain will be kept. For highly autocorrelated
#' chains, that are only examined by few iterations (say less than 1000),
#' the \code{geweke.diag} might fail to detect convergence. In such cases it is
#' essential to look a chain free from autocorelation.
#' @param mn An integer defining the minimum number of individuals per cluster.
#' @param heap A numeric value saying to which value the data might be heaped.
#' Or a list with with such values and names identical to the variables with heaps
#' (see \code{list_of_spikes_maker} for details.)
#' @param rounding_degrees A numeric vector with the presumed rounding degrees. Or a list with rounding degrees,
#' where each list element has the name of a rounded continuous variable. Such a list can be generated
#' using \code{list_of_rounding_degrees_maker(data)}.
#' @param rounding_formula A formula with the model formula for the latent rounding tendency G.
#' Or a list with model formulas for G, where each list element has the name of a rounded continuous variable.
#' Such a list can be generated
#' @return A data.frame where the values, that have a missing value in the original
#' dataset, are imputed.
#' @export
imputationcycle <- function(data_before,
                            original_data,
                            NA_locator,
                            fe, interaction_names,
                            list_of_types,
                            nitt,
                            burnin,
                            thin,
                            mn,
                            heap = 0,
                            rounding_degrees = NULL,
                            rounding_formula){

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

  #update missing rates after the sorting
  missing_rates <- colMeans(NA_locator)

  #get variables with missing values
  incomplete_variables <- names(missing_rates)[missing_rates > 0]

  variables_to_impute <- union(incomplete_variables, names(list_of_types)[list_of_types == "roundedcont" |
                           list_of_types == "interval"])

  #initialize list with the chains of Gibbs-samples
  chains <- list()

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

    #To preserve the relation of y and x properly,
    #y has to be included as covariate in the imputation model for x.
    #But not only as a fix effect covariate,
    #but also as a random effect covariate in the cases where x
    #is a random effects covariate in the analysis model for y.
    if(l2 %in% fe$randomeffects_varname){
      tmp_Z <- cbind(tmp_Z, X[, fe$target_varname])
    }

    tmp_Z <- tmp_Z[, stats::complete.cases(t(tmp_Z)), drop = FALSE]
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

    if(is.list(rounding_degrees)){
      rounding_degrees_tmp <- rounding_degrees[[l2]]
    }else{
      rounding_degrees_tmp <- rounding_degrees
    }

    if(tmp_type == "binary"){
      if(use_single_level){

        #print("Currently the single level model.")
        imp <- imp_binary_single(y_imp = data_before[, l2],
                                 X_imp = tmp_X,
                                 rounding_degrees = rounding_degrees_tmp)

      }else{

        #print("Currently the multilevel model.")
        tmp <- imp_binary_multi(y_imp = data_before[, l2],
                                X_imp = tmp_X,
                                Z_imp = tmp_Z,
                                clID = data_before[, fe$clID_varname],
                                nitt = nitt,
                                burnin = burnin,
                                thin = thin,
                                rounding_degrees = rounding_degrees_tmp)

        imp <- tmp$y_ret
        chains[[l2]]$Sol <- tmp$Sol
        chains[[l2]]$VCV <- tmp$VCV
      }
    }

    if(tmp_type == "cont"){

      if(use_single_level){

        #print("Currently the single level model.")
        imp <- imp_cont_single(y_imp = data_before[, l2],
                               X_imp = tmp_X,
                               rounding_degrees = rounding_degrees_tmp)

      }else{

        #print("Currently the multilevel model.")
        tmp <- imp_cont_multi(y_imp = data_before[, l2],
                              X_imp = tmp_X,
                              Z_imp = tmp_Z,
                              clID = data_before[, fe$clID_varname],
                              nitt = nitt,
                              burnin = burnin,
                              thin = thin,
                              rounding_degrees = rounding_degrees_tmp)
        imp <- tmp$y_ret
        chains[[l2]]$Sol <- tmp$Sol
        chains[[l2]]$VCV <- tmp$VCV
      }

    }

    if(tmp_type == "semicont"){
      if(is.list(heap)){
        heap_tmp <- heap[[l2]]
      }else{
        heap_tmp <- heap
      }

      if(use_single_level){


        #print("Currently the single level model.")
        imp <- imp_semicont_single(y_imp = data_before[, l2],
                                   X_imp = tmp_X,
                                   heap = heap_tmp,
                                   rounding_degrees = rounding_degrees_tmp)
      }else{

        #print("Currently the multilevel model.")
        tmp <- imp_semicont_multi(y_imp = data_before[, l2],
                                  X_imp = tmp_X,
                                  Z_imp = tmp_Z,
                                  clID = data_before[, fe$clID_varname],
                                  heap = heap_tmp,
                                  nitt = nitt,
                                  burnin = burnin,
                                  thin = thin,
                                  rounding_degrees = rounding_degrees_tmp)

        imp <- tmp$y_ret
        chains[[l2]]$Sol <- tmp$Sol
        chains[[l2]]$VCV <- tmp$VCV

      }

    }


    if(tmp_type == "roundedcont"){
      #yet, we don't have a multilevel imputation for rounded incomes

      #Set up the data.frame including the variable explaning G, the latent rounding tendency
      vars_usable <- rounding_formula[[l2]]
      vars_usable <- vars_usable[vars_usable %in% colnames(data_before)]
      PSI <- data_before[, vars_usable, drop = FALSE]

      imp <- imp_roundedcont(y_df = data_before[, l2, drop = FALSE],
                             X = tmp_X,
                             PSI = PSI,
                             rounding_degrees = rounding_degrees_tmp)

    }

    if(tmp_type == "interval"){

      #yet, we don't have a multilevel imputation for interval data
      imp <- imp_interval(y_imp = data_before[, l2],
                          X_imp = tmp_X,
                          rounding_degrees = rounding_degrees_tmp)

    }


    if(tmp_type == "count"){

      if(use_single_level){

        #print("Currently the single level model.")
        tmp <- imp_count_single(y_imp = data_before[, l2],
                                X_imp = tmp_X,
                                nitt = nitt,
                                burnin = burnin,
                                thin = thin,
                                rounding_degrees = rounding_degrees_tmp)

        imp <- tmp$y_ret
        chains[[l2]]$Sol <- tmp$Sol
        chains[[l2]]$VCV <- tmp$VCV


      }else{

        #print("Currently the multilevel model.")
        tmp <- imp_count_multi(y_imp = data_before[, l2],
                               X_imp = tmp_X,
                               Z_imp = tmp_Z,
                               clID = data_before[, fe$clID_varname],
                               nitt = nitt,
                               burnin = burnin,
                               thin = thin,
                               rounding_degrees = rounding_degrees_tmp)

        imp <- tmp$y_ret
        chains[[l2]]$Sol <- tmp$Sol
        chains[[l2]]$VCV <- tmp$VCV

      }

    }

    if(tmp_type == "categorical"){
      data_before[, l2] <- original_data[, l2]
      # if the cluster ID is to be imputed, obviously we cannot run a multilevel model.
      if(use_single_level || l2 == fe$clID_varname){

        #print("Currently the single level model.")
        imp <- imp_cat_single(y_imp = data_before[, l2],
                              X_imp = tmp_X,
                              rounding_degrees = rounding_degrees_tmp)

      }else{

        #print("Currently the multilevel model.")
        tmp <- imp_cat_multi(y_imp = as.factor(data_before[, l2]),
                             X_imp = tmp_X,
                             Z_imp = tmp_Z,
                             clID = data_before[, fe$clID_varname],
                             nitt = nitt,
                             burnin = burnin,
                             thin = thin,
                             rounding_degrees = rounding_degrees_tmp)

        imp <- tmp$y_ret
        chains[[l2]]$Sol <- tmp$Sol
        chains[[l2]]$VCV <- tmp$VCV
      }
    }

    if(tmp_type == "ordered_categorical"){

      if(use_single_level){

        #print("Currently the single level model.")
        imp <- imp_orderedcat_single(y_imp = data_before[, l2],
                                     X_imp = tmp_X,
                                     rounding_degrees = rounding_degrees_tmp)

      }else{

        #print("Currently the multilevel model.")
        tmp <- imp_orderedcat_multi(y_imp = as.factor(data_before[, l2]),
                                    X_imp = tmp_X,
                                    Z_imp = tmp_Z,
                                    clID = data_before[, fe$clID_varname],
                                    nitt = nitt,
                                    burnin = burnin,
                                    thin = thin,
                                    rounding_degrees = rounding_degrees_tmp)

        imp <- tmp$y_ret
        chains[[l2]]$Sol <- tmp$Sol
        chains[[l2]]$VCV <- tmp$VCV
      }
    }

    if(tmp_type == "intercept"){
      imp <- sample_imp(data_before[, l2])
    }
    data_before[, l2] <- imp
  }
  data_after <- data_before[, original_ordering, drop = FALSE]
  ret <- list(data_after = data_after, chains = chains)
  return(ret)
}


#https://stackoverflow.com/questions/5789982/reset-par-to-the-default-values-at-startup
#https://stackoverflow.com/users/429846/gavin-simpson
#' Function to reset all graphics parameters
#'
#' With this function other functions can manipulate the \code{par} settings and
#' when they are finished restore the settings to the state before the function call.
#' The package \code{hmi} uses this mechanism for its function \code{chaincheck}.
#' @author Gavin Simpson
#' @references https://stackoverflow.com/questions/5789982/reset-par-to-the-default-values-at-startup
#' @export
resetPar <- function(){
  grDevices::dev.new()
  op <- graphics::par(no.readonly = TRUE)
  grDevices::dev.off()
  op
}

#' Checking the chains on convergence
#'
#' Formally tests the Gibbs-sampling chains on convergence.
#' After the burn in is discarded, the remaining iterations of each chain
#' are tested following Geweke (1992).
#' In this test, the arithmetic means and their standard errors of the first 10\%
#' and last 50\% of the chain (from now on always after discarding the burn in)
#' are compared. In case of a stationary distribution, both means have the same
#' expected value. The difference between both arithmetic means is divided
#'  THE ONE (CHECK) standard error.
#' This is the Z-score, the test statistic.
#' Chains not passing the test will be plotted.
#' Each plot will flag which (fixed effect or variance) parameter was tested;
#' and what variable was to be imputed and the cycle and imputation run.
#' To see the next plot, the user has to hit <Return> (or "Enter").
#' @param mids A mids object generated by hmi
#' (alternatively a list), having an element called "gibbs" with the chains
#' of the Gibbs-sampler runs.
#' @param alpha A numeric value between 0 and 1 for the desired significance level
#' of the test on convergence.
#' @param plot Logical. Shall the chains be plotted in a traceplot or not.
#' If the number of iterations and cycles is large, click through all traceplots
#' can be interminable.
#' @references J Geweke (1992): Evaluating the accuracy of sampling based approaches
#' to calculating posterior moments. In Bayesian Statistics 4
#' (ed. JB Bernando, JO Berger, AP Dawid and Adrian FM Smith) (pp. 169-193).
#' Clarendon Press, Oxford, UK.
#' @export
chaincheck <- function(mids, alpha = 0.01, plot = TRUE){
  if(!is.numeric(alpha)){
    stop("alpha has to be a numeric value between 0 and 1.")
  }
  if(length(alpha) > 1){
    stop("alpha has to be a single numeric value between 0 and 1 and not a vector.")
  }
  if(alpha < 0 | alpha > 1){
    stop("alpha has to be between 0 and 1.")
  }

  burnin <- mids$gibbs_para[2]
  thin <- mids$gibbs_para[3]

  # by hitting <Return> (a.k.a. [Enter]) the user shall be able to go through the problematic
  # chains by hand
  graphics::par(ask = TRUE)


  counter_all <- 0
  counter_problematic <- 0
  for(i in names(mids$gibbs)){ # the imputations
    for(j in names(mids$gibbs[[i]])){ #the cycles
      for(k in names(mids$gibbs[[i]][[j]])){ #the variables
        for(l in 1: ncol(mids$gibbs[[i]][[j]][[k]]$Sol)){ #the fixed effects parameters
          tmp0 <- as.numeric(mids$gibbs[[i]][[j]][[k]]$Sol[, l])
          tmp1 <- coda::as.mcmc(as.numeric(tmp0))
          tmp2 <- coda::geweke.diag(tmp1, frac1 = 0.1, frac2 = 0.5)$z

          counter_all <- counter_all + 1

          if(abs(tmp2) > stats::qnorm(1 - alpha/2)){
            if(plot){
              coda::traceplot(tmp1, main = paste(i, "; ", j, ";\n",
                                               "Imp. of variable ", k, ";\n fix parameter ", l,
                                               "; z-value: ", round(tmp2, digits = 3),
                                               sep = ""))
            }
            counter_problematic <- counter_problematic + 1

          }
        }

        #In VCV the off-diagonal elements of the random effects covariance matrix are
        #included twice, so those dublicated columns have to be removed.
        vcv <- mids$gibbs[[i]][[j]][[k]]$VCV
        vcv2 <- t(unique(t(as.matrix(vcv))))
        for(l in 1: ncol(vcv2)){ # the variance parameters
          tmp0 <- vcv2[, l]

          tmp1 <- coda::as.mcmc(tmp0)

          tmp2 <- coda::geweke.diag(tmp1, frac1 = 0.1, frac2 = 0.5)$z
          #for the imputation of some variables, some variance parameters have to be fixed
          #(e.g. binary variables). Those would give an NaN for the geweke.diag, so we set the
          #test statistics to 0 in these cases.
          if(get_type(tmp0) == "intercept"){
            tmp2 <- 0
          }
          counter_all <- counter_all + 1

          if(abs(tmp2) > stats::qnorm(1 - alpha/2)){
            if(plot){
              coda::traceplot(tmp1, main = paste(i, "; ", j, ";\n",
                                                 "Imp. of variable ", k, ";\n variance parameter ", l,
                                                 "; z-value: ", round(tmp2, digits = 3),
                                                 sep = ""))
            }
            counter_problematic <- counter_problematic + 1

          }

        }
      }
    }
  }
  cat(paste(counter_problematic, " out of ", counter_all,
            " chains (",round(counter_problematic/counter_all * 100, digits = 2),"%) " ,
            "did not pass the convergence test.\n",
            "For alpha = ", alpha, ", the expected number is ", counter_all*alpha, ".", sep = ""))
  graphics::par(resetPar())
}
