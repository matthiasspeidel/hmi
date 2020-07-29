#' Sample imputation.
#'
#' Function to sample values in a variable from other (observed) values in this variable.
#' So this imputation does not use further covariates.
#' @param variable A vector of size \code{n} with missing values.
#' @return A n x 1 data.frame with the observed and imputed data
#' @examples
#' set.seed(123)
#' sample_imp(c(1, NA, 3, NA, 5))
#' @export
sample_imp <- function(variable){

  if(is.data.frame(variable)){
    stop("You passed a data.frame instead of a vector to sample_imp.")
  }
  if(all(is.na(variable))) stop("Variable consists only of NAs.")

  ret <- data.frame(target = variable, stringsAsFactors = FALSE)
  need_replacement <- is.na(variable) | is.infinite(variable)
  ret[need_replacement, 1] <- sample(size = sum(need_replacement),
                             variable[!need_replacement], replace = TRUE)
  return(ret)
}


#' Get the type of variables.
#'
#' Function checks of which type a variable is. The types are listed below
#' (together with a rough summary of the classification mechanism).
#' \itemize{
#'  \item continuous (numeric values, or integers with more than 20 different values),
#'  \item semicontinuous (numeric values with more than 10\% of them share the same value),
#'  \item rounded continuous (if more than 50\% of the observations of a continuous variable
#'   are divisible by some rounding degrees)
#'  \item count data (integers).
#'  \item an intercept (the same value for all observations),
#'  \item binary (two different values - like 0s and 1s or "m" and "f"),
#'  \item categorical (the variable is a factor or has more than 3 different values)
#'  \item ordered categorical (the categorical variable is ordered.)
#'}
#'
#' @param variable A variable (vector) from your data set.
#' @param spike A numeric value, denoting the presumed spike of semi-continuous variables.
#' @param rounding_degrees A numeric vector with the presumed rounding degrees.
#' If the rounding_degrees are set to be \code{NULL}, the proposed rounding degrees from the function
#' \code{suggest_rounding_degrees} are used.
#' @return A character denoting the type of \code{variable}.
#' @examples get_type(iris$Sepal.Length); get_type(iris$Species)
#' @export
get_type <- function(variable, spike = NULL,
                     rounding_degrees = c(1, 10, 100, 1000)){

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

    #If too many categories are available, the variable is considered to be continuous.
    if(length(unique(variable)) > 20){
      type <- "cont"
    }
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

      # if there is no difference between a variable and their rounding result,
      # it is considered to be an integer...
      if(max(abs(variable - round(variable)), na.rm = TRUE) == 0){
        type <- "count"

        #... but not if too many categories are available...
        if(length(table(variable)) > 20){
          type <- "cont"
        }

        #... or it is a rounded continuous variable:
        if(sum(apply(outer(decompose_interval(interval = variable)[, "precise"], setdiff(rounding_degrees, 1), '%%') == 0,
                     1, any), na.rm = TRUE)/
           length(variable[!is.na(variable)]) > 0.5){
          type <- "roundedcont"
        }

        return(type)
      }

      #if the variable is numeric and more than 10 \% of the variable values share the same value,
      #We consider this variable as semi-continious.

      if(is.null(spike)){
        number_at_spike <- max(table(variable))
      }else{
        number_at_spike <- sum(variable == spike, na.rm = TRUE)
      }
      if(number_at_spike/length(variable[!is.na(variable)]) > 0.1){
        type <- "semicont"
      }

      # if more than 50 \% of the data are divisible by on of the given rounding degrees,
      # they are considered to be rounded continuous
      # Observed 0s shall not be counted as rounded,
      #because otherwise semi-continuous variables might be considered to be rounded-continuous.

      if((sum(apply(outer(decompose_interval(interval = variable)[, "precise"], setdiff(rounding_degrees, 1), '%%') == 0,
                    1, any), na.rm = TRUE) - sum(variable == 0, na.rm = TRUE))/
        length(variable[!is.na(variable)]) > 0.5){
        type <- "roundedcont"
      }


      return(type)

    }

    if(type == "unknown"){
      type <- "categorical"
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
      #if more than 50 \% of the precise part of an interval variable is rounded,
      #the whole variable is considered to be a rounded continous variable
      tmp <- decompose_interval(interval = variable)[, "precise"]
      if(sum(!is.na(tmp)) == 0) return(type)
      if(sum(apply(outer(tmp, setdiff(rounding_degrees, 1), '%%') == 0,
                   1, any), na.rm = TRUE)/
         length(tmp[!is.na(tmp)]) > 0.5){
        type <- "roundedcont"
      }
      return(type)
    }

    if(ncol(variable) == 1){
      ret <- get_type(variable[, 1], spike = spike, rounding_degrees = rounding_degrees)
      return(ret)
    }

    if(ncol(variable) > 1){
      if(is.matrix(variable)){
        variable <- as.data.frame(variable)
      }

      ret <- array(dim = ncol(variable))

      for(i in 1:length(ret)){
        if(is.list(rounding_degrees)){
          rounding_degrees_tmp <- rounding_degrees[[colnames(variable)[i]]]
        }else{
          rounding_degrees_tmp <- rounding_degrees
        }

        if(is.list(spike)){
          spike_tmp <- spike[[colnames(variable)[i]]]
        }else{
          spike_tmp <- spike
        }

        ret[i] <- get_type(variable[, i], spike = spike_tmp, rounding_degrees = rounding_degrees_tmp)
      }

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

#' Helps the user to make a list of spikes.
#'
#' In \code{hmi} the user can add a list of spikes. This function gives her/him a framework
#' with suggestions. Of course the user can make changes by herself/himself afterwards.
#' For example, the function might wrongly classify a variable to have a spike.
#' @param data the data.frame also passed to \code{hmi}.
#' @return a list with suggested spikes. Each list element has the name of a spiked variable
#' in the data.frame. The elements contain a single numeric denoting the spike found for that variable.
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
#' For example, if a continuous variable has only two observations left, then get_type
#' interpret this as a binary variable and not a continuous.
#' @param data the data.frame also passed to \code{hmi}.
#' @param spike A numeric value or list saying which value in the semi-continuous data might be the spike.
#' @param rounding_degrees A numeric vector with the presumed rounding degrees.
#' @return a list with suggested types. Each list element has the name of a variable
#' in the data.frame. The elements contain a single character denoting the type of the variable.
#' See \code{get_type} for details about the variable types.
#' @export
list_of_types_maker <- function(data, spike = NULL, rounding_degrees = NULL){

  if(!is.data.frame(data)) stop("We need the data in the data.frame format!")
  if(nrow(data) == 0){
    ret <- list()
    for(l2 in colnames(data)){
      ret[[l2]] <- "unknown"
    }
  }
  if(ncol(data) == 0){
    return(NULL)
  }
  if(is.null(data)) return(NULL)
  rounding_degrees_tmp <- rounding_degrees

  ret <- list()

  for(l2 in colnames(data)){
    if(is.list(spike)){
      spike_tmp <- spike[[l2]]
    }else{
      spike_tmp <- spike
    }

    if(is.list(rounding_degrees)){
      rounding_degrees_tmp <- rounding_degrees[[l2]]
    }

    #The specifications of variables depends on 4 siutations.
    #Setting 1: rounding_degrees is NULL
    #Setting 2: rounding_degrees is a vector
    #Setting 3a: rounding_degrees is list, with no element for the current variable
    #Setting 3b: rounding_degrees is list, with a specific vector for the current variable

    #The consequences for these settings are:
    #Setting 1: Test using rounding degrees 1, 10, 100, 1000
    #Setting 2: Test using the vector
    #Setting 3a: Test using rounding degrees 1, 10, 100, 1000
    #Setting 3b: define as roundedcont
    if(is.null(rounding_degrees)){       ## Setting 1 ##
      rounding_degrees_tmp <- c(1, 10, 100, 1000)
      tmp_type <- get_type(data[, l2],
                           spike = spike_tmp, rounding_degrees = rounding_degrees_tmp)

    }else if(is.integer(rounding_degrees) | is.numeric(rounding_degrees)){      ## Setting 2 ##
      rounding_degrees_tmp <- rounding_degrees
      tmp_type <- get_type(data[, l2],
                           spike = spike_tmp, rounding_degrees = rounding_degrees_tmp)

    }else if(is.list(rounding_degrees) & is.null(rounding_degrees[[l2]])){      ## Setting 3a ##
      rounding_degrees_tmp <- c(1, 10, 100, 1000)
      tmp_type <- get_type(data[, l2],
                           spike = spike_tmp, rounding_degrees = rounding_degrees_tmp)
    }else if(is.list(rounding_degrees) & !is.null(rounding_degrees[[l2]])){      ## Setting 3b ##
      tmp_type <- "roundedcont"
      rounding_degrees_tmp <- rounding_degrees[[l2]]
    }

    #If the type is still NULL, set it to sample_imp
    if(is.null(tmp_type)){
      tmp_type <- "sample_imp"

    }

    ret[[l2]] <- tmp_type
  }

  return(ret)
}

#' Function to get all factors
#'
#' Function to get all factors (not limited to prime factors) of an integer.
#' @param x A single integer; no vector.
#' @return A numeric vector with the factors
#' @references based on stackoverflow.com/questions/6424856 "R Function for returning ALL factors"
#'  answer by Chase
factors <- function(x){
  if(!is.numeric(x)) return(NA)
  if(length(x) != 1) return(NA)
  if(is.na(x)) return(NA)
  if(is.infinite(x)) return(NA)
  if(x %% 1 != 0) return(NA)
  x <- as.integer(x)
  div <- seq_len(abs(x))
  return(div[x %% div == 0L])
}


#' suggesting rounding degrees
#'
#' A function that suggests some rounding degrees of a continuous variable
#' (classically formatted or as interval object).
#' The basic idea is 1. to count which factor is observed in the data more often than expected.
#' 2. to check whether a factor can explain at least three observed piles of observations in the data
#' 3. to check whether a factor explains at least 20 \% of observations (additional to previous factors).
#' Factors fulfilling this premises are returned as suggested rounding degrees.
#' @param x A vector or \code{interval} object.
suggest_rounding_degrees <- function(x){

  if(!(is.vector(x) | is_interval(x))){
    return(NULL)
  }

  if(is.factor(x)){
    return(NULL)
  }
  variable <- decompose_interval(interval = x)[, "precise"]
  tmpvariable <- variable
  n <- length(variable)
  tab_factors <- table(unlist(sapply(tmpvariable, factors)))/n

  if(length(tab_factors) == 0){
    return(NULL)
  }else{

    values <- as.numeric(names(tab_factors))
    observed_values <- setdiff(as.numeric(names(table(variable))), 0)
    #for later, only those values shall be considered if the observed number of indiduvials,
    #rounded to this degree, exceeds the expected number at least by two.
    #Example: from 10000 individuals, 2000 are expected to round to 5.
    #If 4000 or more indivudials are observed 5 is considered to be a possible rounding degree
    candidates_values <- values[values * tab_factors >= 2]
    #instead of (absolute freq) / (n/s) > 2, it can be used
    #(absolute freq) * s/n > 2  which is equivalent to
    #(relative freq) * s > 2.
    if(length(candidates_values) == 0) return(NULL)

    ret <- array(dim = 0)
    counter <- 1
    for(cv in sort(candidates_values, decreasing = TRUE)){
      if(sum(observed_values %% cv == 0) >= 3){

        #check which portion of heaping could be explained by cv.
        tmp <- tab_factors[names(tab_factors) == cv]
        if(length(tmp) > 0 && tmp >= 0.2){
          ret[counter] <- cv
          counter <- counter + 1

          #update tab_factors: every following (smaller) factor which can divide cv without rest
          #has to explain 20 \% more of the data than cv. For this reason, observation divisable by
          #cv are removed from the tmpvariable
          tmpvariable <- tmpvariable[tmpvariable %% cv != 0]
          #update tab_factors
          tab_factors <- table(unlist(sapply(tmpvariable, factors)))/n

        }
      }
    }
    suggested_rounding_degrees <- sort(unique(c(1, ret)))
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


#' Helps the user to make a list of rounding formulas for the rounding degrees
#'
#' In \code{hmi} the user can add a list of rounding formulas for each variable suffering rounding.
#' This function gives him/her a framework with suggestions. Of course the user can make changes
#' by herself/himself afterwards.
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
#' This function applies the analysis the user wants to run on every imputed dataset.
#' The results from every dataset are pooled by simply averaging. So the user has to make sure that
#' averaging the analysis produces results is meaningful. Currently variance estimates for the averaged
#' results are not implemented.
#' @param mids A \code{mids} (multiply imputed data set) object.
#' Either from the \code{hmi} imputation function or \code{mice}.
#' @param analysis_function A user generated function that gets a completed data set,
#' runs the model and returns all model parameters
#' he or she is interested in in a vector. See examples below.
#' @return A vector with all averaged results.
#' @examples
#' \dontrun{
#'  data(Gcsemv, package = "hmi")
#'
#'  model_formula <- written ~ 1 + gender + coursework + (1 + gender|school)
#'
#'  set.seed(123)
#'  dat_imputed <- hmi(data = Gcsemv, model_formula = model_formula, M = 2, maxit = 2)
#'
#'  my_analysis <- function(complete_data){
#'   # In this list, you can write all the parameters you are interested in.
#'   # Those will be averaged.
#'   # So make sure that averaging makes sense and that you only put in single numeric values.
#'   parameters_of_interest <- list()
#'
#'   # ---- write in the following lines, what you are interetest in to do with your complete_data
#'   # the following lines are an example where the analyst is interested in the fixed intercept
#'   # and fixed slope and the random intercepts variance,
#'   # the random slopes variance and their covariance
#'   my_model <- lmer(model_formula, data = complete_data)
#'
#'   parameters_of_interest[[1]] <- fixef(my_model)
#'   parameters_of_interest[[2]] <- lme4::VarCorr(my_model)[[1]][,]
#'   ret <- unlist(parameters_of_interest)# This line is essential if the elements of interest
#'   #should be labeled in the following line.
#'   names(ret) <-
#'     c("beta_intercept", "beta_gender", "beta_coursework", "sigma0", "sigma01", "sigma10", "sigma1")
#'
#'   return(ret)
#' }
#' hmi_pool(mids = dat_imputed, analysis_function = my_analysis)
#' #if you are interested in fixed effects only, consider pool from mice:
#' pool(with(data = dat_imputed, expr = lmer(written ~ 1 + gender + coursework + (1 + gender|school))))
#' }
#' @importFrom mice complete
#' @export
hmi_pool <- function(mids, analysis_function){

  if (!mice::is.mids(mids)){
    stop("Your multiple imputed data set must have class mids.")
  }

  results <- list()

  for (i in 1:mids$m) {
    results[[i]] <- unlist(analysis_function(mice::complete(mids, i)))
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
#' @param parnames A character vector with the names of the elements in para.
#' @param X_in_negloglik The data.frame of covariates explaining Y, the observed target variable.
#' It has to has n rows (with n being the number of precise, imprecise and missing observations).
#' @param PSI_in_negloglik The data.frame of covariates explaining G, the latent rounding tendency.
#' Without the target variable.
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
#' @param rounding_degrees A numeric vector with the presumed rounding degrees for Y.
#' @references Joerg Drechsler, Hans Kiesl, Matthias Speidel (2015):
#' "MI Double Feature: Multiple Imputation to Address Nonresponse and Rounding Errors in Income Questions",
#' Austrian Journal of Statistics, Vol. 44, No. 2, \url{http://dx.doi.org/10.17713/ajs.v44i2.77}
#' @return An integer equal to the (sum of the) negative log-likelihood contributions (of the observations)
negloglik <- function(para,
                      parnames = names(para),
                      X_in_negloglik,
                      PSI_in_negloglik,
                      y_precise_stand,
                      lower_bounds = NA, upper_bounds = NA,
                      my_g, sd_of_y_precise,
                      indicator_precise,
                      indicator_imprecise,
                      indicator_outliers,
                      rounding_degrees = c(1, 10, 100, 1000)){

  names(para) <- parnames

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
    reg_coefs <- matrix(para[grep("^coef_y_on_x", names(para))], ncol = 1)
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
    as.matrix(PSI_in_negloglik[!indicator_outliers, , drop = FALSE]) %*% G_coefs

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

  #Sum up the likelihood contributions according to the rounding degree:
  #For observations rounded to the first rounding degree only the first likelihoodcontribution is considered
  #For observations rounded to the fourth rounding degree, the first four contributions are added up.
  # Observations with no rounding, are excluded.
  result_try <- apply(merged, 1, function(x) sum(x[1:x[length(x)]]))
  result_try[is.na(result_try)] <- 0
  result_try <- result_try[my_g[!indicator_outliers] != 0]
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
#' @param parnames A character vector with the names of the elements in para.
#' @param X the data.frame of covariates.
#' @param lower_bounds the lower bound of an interval variable.
#' @param upper_bounds the upper bound of an interval variable.
#' @references Joerg Drechsler, Hans Kiesl, Matthias Speidel (2015):
#' "MI Double Feature: Multiple Imputation to Address Nonresponse and Rounding Errors in Income Questions",
#' Austrian Journal of Statistics, Vol. 44, No. 2, \url{http://dx.doi.org/10.17713/ajs.v44i2.77}
#' @return An integer equal to the (sum of the) negative log-likelihood contributions (of the observations)
negloglik2_intervalsonly <- function(para, parnames = names(para), X, lower_bounds, upper_bounds){

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
                             variable_names_in_data = colnames(data), data){


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

  #If a dot was used in the model_formula (e.g. y ~ 0 + .) stats::terms() need a data.frame as all variables from this dataset will be used.
  #In this case, here an artificial dataset with a (hopefully) unique variable name is used.
  fixedeffects_varname <- all.vars(stats::delete.response(
    stats::terms(lme4::nobars(model_formula), data = data.frame("ThisVariableWasIntroducedByHMIBecauseADotAppearedInTheFormula" = 1:2))
    ))
  # --check if intercept is present
  fixed_intercept_exists <- attributes(stats::delete.response(
    stats::terms(lme4::nobars(model_formula), data = data.frame("ThisVariableWasIntroducedByHMIBecauseADotAppearedInTheFormula" = 1:2))
    ))$intercept == 1

 # if(fixedeffects_varname == "ThisVariableWasIntroducedByHMIBecauseADotAppearedInTheFormula"){
#    fixed_intercept_exists <- FALSE
#  }
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
  #if there are no random effects at all, there definitely is random intercept
  if(length(randomeffects_varname) == 0) return(FALSE)
  #if the only random effect variable is "0", the input makes no sense, and of course no random intercept is specified
  if(randomeffects_varname == "0") return(FALSE)
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
#' @param sna Boolean: if \code{TRUE}, \code{NA}s are kept as standard \code{NA}s.
#' Otherwise they are turned into \code{"-Inf;Inf"}.
#' @seealso \link[hmi]{generate_interval}
#' @return A vector of class \code{interval}.
#' @examples as.interval(c("1000;2000", "700;700", NA))
#' @export
as.interval <- function(x, sna = FALSE){
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
  if(sna) orgNA <- is.na(x)

  x[is.na(x)] <- c("NA;NA")
  raw <- unlist(strsplit(x, ";"))
  raw[raw == "NA"] <- NA
  lower <- as.numeric(raw[seq(1, 2*length(x), by = 2)])
  upper <- as.numeric(raw[seq(2, 2*length(x), by = 2)])

  if(sna){
    lower[orgNA] <- NA
    upper[orgNA] <- NA
    ret <- generate_interval(lower, upper, sna = TRUE)
  }else{
    ret <- generate_interval(lower, upper)
  }

  return(ret)
}

#' Function to generate an interval object
#'
#' Function to generate transform into an
#' \code{interval} object from numeric (or character) vectors.
#' @param lower a vector with the lower values of a variable
#' @param upper a vector with the upper values of a variable
#' @param sna Boolean: if \code{TRUE}, \code{NA}s are kept as standard \code{NA}s.
#' Otherwise they are turned into \code{"-Inf;Inf"}.
#' @return a character vector with the lower and upper bound values.
#' @export
generate_interval <- function(lower, upper, sna = FALSE){
  if(length(lower) != length(upper)) stop("lower and upper need the same length")

  orgNAlower <- is.na(lower)
  orgNAupper <- is.na(upper)
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
  if(sna){
    ret[orgNAlower & orgNAupper] <- NA
  }
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
  reg_low <- "^((-){0,1}[[:digit:]]{1,}(\\.[[:digit:]]{1,}){0,1}|(-){0,1}Inf|NA)"

  #the end of the interval is one or more digits optionally preceded by an "-"
  # optionally followed by one period (or "full stop") and one or more digits
  #or the interval ends with "Inf"
  reg_up <- "((-){0,1}[[:digit:]]{1,}(\\.[[:digit:]]{1,}){0,1}|(-){0,1}Inf|NA)$"
  #currently no check on whether the interval is valid (lower bound does not exceed upper bound) is made.
  matches <- grep(paste(reg_low, reg_up, sep = ";"), x)

  #currently c("100", "100;200") is not classified as interval. This might be changed in future by setting
  #return(length(matches) >= 1)
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
#' @export
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

  interval <- as.character(interval)

  interval[is.na(interval)] <- "NA;NA"
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

  if(!(is.vector(interval) | is.factor(interval) | is.numeric(interval)| is_interval(interval))){
    stop("interval has to be a vector, factor, numeric or interval.")
  }

  if(is.factor(interval)){
    interval <- as.character(interval)
  }

  if(is.character(interval)){
    interval <- as.interval(interval)
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

  # Get the precise observations of the interval (the imprecise observations will be NA)
  precise <- ifelse(abs(dists) < 1e-20, tmp[, 1], NA)

  # Get the imprecise observations of the interval (the precise will be NA)
  imprecise <- as.interval(ifelse(abs(dists) >= 1e-20 | is.na(dists), interval, NA), sna = TRUE)

  # Split the imprecise observations up into the lower and upper bound
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
#' This function replaces observations with "-Inf;Inf" in an interval,
#' which basically means "no information available", with the standard NAs
#' (therefore the name 'sna').
#' Observations with a finite bound (e.g.x = "0;Inf") are not replaced,
#' since they contain information (here: "x is positive").
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
#' @return A \code{data.frame} where the interval variables are stored as \code{interval} objects.
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
#' @param sort A character specifying how the values should be sorted if only one variable is to be plotted.
#' By default they are sorted according to their position in the data set.
#' \code{sort = "lowerbound_increasing"} sorts the data primarily by their lower bound, and secondarily
#' (this means for equal lower bounds) by their upper bounds. Both in increasing order.
#' For \code{sort = "lowerbound_decreasing"}, both happens in decreasing order.
#' \code{sort = "mostprecise_increasing"} sorts the data by their length of the interval they represent,
#' and within equal lengths by the lower bound. Both in increasing order.
#' For \code{sort = "mostprecise_decreasing"}, both happens in decreasing order.
#' @param ... graphical parameters such as \code{main}.
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
    if(sort == "lowerbound_increasing"){
      tmp <- decompose_interval(y)
      y_split <- y_split[order(tmp[, "lower_general"], tmp[, "upper_general"]),]
      if(is.null(xlab)) myxlab <- "Index of sorted values"
    }

    if(sort == "lowerbound_decreasing"){
      tmp <- decompose_interval(y)
      y_split <- y_split[order(tmp[, "lower_general"], tmp[, "upper_general"], decreasing = TRUE),]
      if(is.null(xlab)) myxlab <- "Index of sorted values"
    }

    if(sort == "mostprecise_increasing"){
      interval_lengths <- decompose_interval(y)[, "upper_general"] -
        decompose_interval(y)[, "lower_general"]

      y_split <- y_split[order(interval_lengths, y_split[, 1]),]
      if(is.null(xlab)) myxlab <- "Index of sorted values"
    }

    if(sort == "mostprecise_decreasing"){
      interval_lengths <- decompose_interval(y)[, "upper_general"] -
        decompose_interval(y)[, "lower_general"]

      y_split <- y_split[order(interval_lengths, y_split[, 1], decreasing = TRUE),]
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

#' Tabulating interval objects
#'
#' Function to tabulate interval objects
#' @param x In its most save way, \code{x} is an object from class \code{interval}.
#' @param ... Other parameters passed to \code{table}.
#' @return A table.
#' @export
table <- function(x, ...) {
  UseMethod('table', x)
}

#' @param sort A character specifying how the values should be sorted if only one variable is to be plotted.
#' \code{sort = "lowerbound_increasing"} (the default) sorts the data primarily by their lower bound, and secondarily
#' (this means for equal lower bounds) by their upper bounds. Both in increasing order.
#' For \code{sort = "lowerbound_decreasing"} both happens in decreasing order.
#' \code{sort = "mostprecise_increasing"} sorts the data by their length of the interval they represent,
#' and within equal lengths by the lower bound. Both in increasing order.
#' \code{sort = "mostprecise_decreasing"} both happens in decreasing order.
#' @rdname table
#' @export
table.interval <- function(x, sort = "lowerbound_increasing",  ...){
  tab <- table(unclass(x), ...)
  tmp <- decompose_interval(as.interval(names(tab)))

  if(sort == "lowerbound_increasing"){
    tab <- tab[order(tmp[, "lower_general"], tmp[, "upper_general"])]
  }

  if(sort == "lowerbound_decreasing"){
    tab <- tab[order(tmp[, "lower_general"], tmp[, "upper_general"], decreasing = TRUE)]
  }

  if(sort == "mostprecise_increasing"){
    interval_lengths <- tmp[, "upper_general"] -
      tmp[, "lower_general"]
    tab <- tab[order(interval_lengths, tmp[, "lower_general"])]
  }

  if(sort == "mostprecise_decreasing"){
    interval_lengths <- tmp[, "upper_general"] -
      tmp[, "lower_general"]
    tab <- tab[order(interval_lengths, tmp[, "lower_general"], decreasing = TRUE)]
  }

  return(tab)
}

#' @rdname table
#' @export
table.default <- function(x, ...) {
  return(base::table(x, ...))
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
  return(generate_interval(lower = tmp[, 1], upper = tmp[, 2]))
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

#' Modulo function
#'
#' Modulo function for interval objects
#' @param interval an object from class interval
#' @param x an single element or vector to add to the interval object
#' @return a vector with the modulo for the precise elements in \code{interval}.
#' For imprecise elements, \code{NA} is returned.
#' @export
#' @rdname interval-modulo
#' \alias{\%<unescaped bksl>\%}
"%%.interval" <- function(x, interval){
  return(decompose_interval(interval)[, "precise"] %% x)
}


#' Log function for interval objects
#'
#' Log function for interval objects
#' @param x numeric vector or interval object
#' @param ... further arguments passed to \code{log}
#' @export
log.interval <- function(x, ...){
  tmp <- base::log(split_interval(x, ...))
  return(generate_interval(lower = tmp[, 1], upper = tmp[, 2]))
}

#' Exp function for interval objects
#'
#' Exp function for interval objects
#' @param x numeric vector or interval object
#' @param ... further arguments passed to \code{exp}.
#' @export
exp.interval <- function(x, ...){
  tmp <- base::exp(split_interval(x, ...))
  return(generate_interval(lower = tmp[, 1], upper = tmp[, 2]))
}


#' Power function
#'
#' Taking the power of interval objects
#' @param interval an object from class interval
#' @param x an single numeric to potentiate
#' @return an interval object
#' @export
#' @rdname interval-power
"^.interval" <- function(interval, x){
  tmp <- split_interval(interval) ^ x
  generate_interval(lower = tmp[, 1], upper = tmp[, 2])
}


#' Sqrt function for interval objects
#'
#' Sqrt function for interval objects
#' @param x numeric vector or interval object
#' @param ... further arguments passed to \code{sqrt}.
#' @export
sqrt.interval <- function(x, ...){
  tmp <- base::sqrt(split_interval(x, ...))
  return(generate_interval(lower = tmp[, 1], upper = tmp[, 2]))
}


#' Round function for interval objects
#'
#' Round function for interval objects
#' @param x numeric vector or interval object
#' @param ... further arguments passed to \code{round}.
#' @export
round.interval <- function(x, ...){
  tmp <- base::round(split_interval(x, ...))
  return(generate_interval(lower = tmp[, 1], upper = tmp[, 2]))
}

#' Floor function for interval objects
#' Floor function for interval objects
#' @param x numeric vector or interval object
#' @param ... further arguments passed to \code{floor}.
#' @export
floor.interval <- function(x, ...){
  tmp <- base::floor(split_interval(x, ...))
  return(generate_interval(lower = tmp[, 1], upper = tmp[, 2]))
}


#' Ceiling funtion for intervals
#'
#' Ceiling funtion for intervals
#' @param x numeric vector or interval object
#' @param ... further arguments passed to \code{ceiling}.
#' @export
ceiling.interval <- function(x, ...){
  tmp <- base::ceiling(split_interval(x, ...))
  ret = generate_interval(lower = tmp[, 1], upper = tmp[, 2])
  return(ret)
}


#' Head for intervals
#'
#' Head function for intervals returning the first elements of an \code{interval} object
#' @param x vector, matrix, table, data.frame or interval object
#' @param ... further arguments passed to \code{head}.
head.interval <- function(x, ...){
  tmp <- utils::head(split_interval(x, ...))
  return(generate_interval(lower = tmp[, 1], upper = tmp[, 2]))
}

#' Tail for intervals
#'
#' Tail function for intervals returning the last elements of an \code{interval} object
#' @param x vector, matrix, table, data.frame or interval object
#' @param ... further arguments passed to \code{tail}.
tail.interval <- function(x, ...){
  tmp <- utils::tail(split_interval(x, ...))
  return(generate_interval(lower = tmp[, 1], upper = tmp[, 2]))
}

#' Index for interval
#'
#' Function to index elements of an interval object
#' @param obj the interval object
#' @param index the index of the elements to replace
`[.interval`  <- function(obj, index){
  ret <- as.interval(as.character(obj)[index], sna = TRUE)
  return(ret)
}


#' Replace for interval
#'
#' Function to replace elements of an interval object
#' @param obj the interval object
#' @param index the index of the elements to replace
#' @param value the value the replaced elements shall take
`[<-.interval` <- function(obj, index, value){
  ret <- as.character(obj)
  ret[index] <- value
  return(as.interval(ret, sna = TRUE))
}

#' Function to give the center of the interval
#'
#'  Function to give the center of the interval object
#' @param interval an object from class interval
#' @param inf2NA logical. If \code{TRUE}, entries containing -Inf or Inf, will return NA.
#' @return A numeric vector
#' @export
center_interval <- function(interval, inf2NA = FALSE){
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
  if(ncol(X) == 0) return(X)
  types <- array(NA, dim = ncol(X))
  for(i in 1:length(types)){
    types[i] <- get_type(X[, i])
  }
  need_stand_X <- types %in% c("cont", "count", "roundedcont", "semicont")
  X_stand <- X
  tmp <- as.data.frame(scale(X[, need_stand_X]))
  X_stand[, need_stand_X] <- matrix(tmp, ncol = ncol(tmp)) # this avoids having attributes delivered by scale().
  return(as.data.frame(X_stand))
}

#' cleanup data.frames
#'
#' Function to exclude (factor or interval) variables that have too many levels
#' (as they may cause numerical problems), to change binary factors to 0-1 coding
#' (as such factors might generate linear dependent variables) and
#' to remove multiple intercepts.
#' @param X A n times p data.frame with p fixed (or random) effects variables.
#' @param k An integer defining the allowed maximum of levels in a factor covariate.
#' @return A n times (p-r) data.frame, with r being the number of variables with too many factors.
#' @export
cleanup <- function(X, k = Inf){

  if(!is.data.frame(X)) stop("X has to be a data.frame.")

  if(!is.numeric(k)) stop("k has to be a single numeric value.")
  if(length(k) > 1)  stop("k has to be a single numeric value.")

  # work out the type of the variables
  types <- array(dim = ncol(X))
  for(j in 1:length(types)) types[j] <- get_type(X[, j])

  # currently, R interprets interval data as categorical.
  # This can be unfeasible as there could be as many categories as observations.
  # Three practical approaches are 1. removing interval data. 2. splitting them up into
  # two separate variables. 3. Interprete them as factors. Currently we chose option 3.
  # The code in comments can be activated, when approach 2 shall be used.
  interval <- types == "interval"
  #added_variables <- as.data.frame(matrix(NA, nrow = nrow(X), ncol = 2 * sum(interval)))
  #interval_counter <- 1
  for(l in which(interval)){
    X[, l] <- as.factor(X[, l])
    #lower <- sample_imp(split_interval(X[, l])[, 1])
    #upper <- sample_imp(split_interval(X[, l])[, 2])
    #added_variables[, 2 * (interval_counter - 1) + 1:2] <- cbind(lower, upper)

    #interval_counter <- interval_counter + 1
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
  categorical <- types == "categorical" | types == "interval"

  # get the name of the variables, that have too many levels
  too_many_levels <- colnames(X[, categorical, drop = FALSE])[
    apply(X[, categorical, drop = FALSE], 2, function(x) nlevels(factor(x))) > k]

  more_than_10_levels <- colnames(X[, categorical, drop = FALSE])[
    apply(X[, categorical, drop = FALSE], 2, function(x) nlevels(factor(x))) > 10]

  #if(length(more_than_10_levels) > 0){
  #  warning(paste("More than 10 categories found in ", paste(more_than_10_levels, collapse = ", "),
  #                ".\n To prevent instable results, consider setting alpha = 0.2.
  #                \n See ?hmi for a description of alpha.", sep = ""))
  #}


  tmp <- X[, !names(X) %in% too_many_levels, drop = FALSE] # If interval variables shall be split up
  #add  & !interval

  # When multiple intercept variables are present, keep only the first
  intercepts <- types == "intercept"
  if(sum(intercepts) >= 2){
    #in other words: remove all intercepts variables, except the first on
    tmp <- tmp[, -which(intercepts)[-1], drop = FALSE]
  }

  ret <- tmp #currentlich not in use: cbind(tmp, added_variables)
  if(ncol(ret) == 0) stop("No variables left after removing excessive variables.")

  return(ret)
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
#' the standard error.
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
#' @param thin An integer to set the thinning interval range. If thin = 1,
#' every iteration of the Gibbs-sampling chain will be kept. For highly autocorrelated
#' chains, that are only examined by few iterations (say less than 1000),
#' the \code{geweke.diag} might fail to detect convergence. In such cases it is
#' essential to look at a chain free from autocorrelation. When setting \code{thin = NULL},
#' the function will use internally a thinning of \code{max(1, round((nitt-burnin)/1000))}
#' to get approximately 1000 iterations to be tested.
#' @param plot Logical. Shall the chains be plotted in a traceplot or not.
#' If the number of iterations and cycles is large, click through all traceplots
#' can be interminable.
#' @references J Geweke (1992): Evaluating the accuracy of sampling based approaches
#' to calculating posterior moments. In Bayesian Statistics 4
#' (ed. JB Bernando, JO Berger, AP Dawid and Adrian FM Smith) (pp. 169-193).
#' Clarendon Press, Oxford, UK.
#' @export
chaincheck <- function(mids, alpha = 0.01, thin = 1, plot){

  if(rlang::is_missing(plot)){
    if(interactive()){
      plot <- TRUE
    }else{
      plot <- FALSE
    }
  }



  if(!is.numeric(alpha)){
    stop("alpha has to be a numeric value between 0 and 1.")
  }
  if(length(alpha) > 1){
    stop("alpha has to be a single numeric value between 0 and 1 and not a vector.")
  }
  if(alpha < 0 | alpha > 1){
    stop("alpha has to be between 0 and 1.")
  }

  nitt <- mids$gibbs_para[1]
  burnin <- mids$gibbs_para[2]
  #old: thin <- mids$gibbs_para[3]

  #If no thinning parameter was given, it is set to a value that about 1000 iterations are to be
  #examinded
  if(is.null(thin)){
    thin <- max(1, round((nitt-burnin)/1000))
    # by using max(1, ...) we make sure that for example when nitt = 600 and burnin = 200, the thinning is sensible.
  }

  #If no MCMCglmm method was used to impute the data, the lists in mids$gibbs are empty.
  #This is checked by this call.
  if(all(unlist(lapply(mids$gibbs, function(x) all(unlist(lapply(x,  function(y) length(y) == 0 ))))))){
    return("No MCMC chains were found. Probably because only single level imputation methods were used.")
  }

  # by hitting <Return> (a.k.a. [Enter]) the user shall be able to go through the problematic
  # chains by hand
  opar <- graphics::par('ask')
  on.exit(graphics::par(ask = opar))

  graphics::par(ask = TRUE)

  counter_all <- 0
  counter_problematic <- 0
  for(i in names(mids$gibbs)){ # the imputations
    for(j in names(mids$gibbs[[i]])){ #the cycles
      for(k in names(mids$gibbs[[i]][[j]])){ #the variables
        for(l in 1: ncol(mids$gibbs[[i]][[j]][[k]]$Sol)){ #the fixed effects parameters
          tmp0 <- as.numeric(mids$gibbs[[i]][[j]][[k]]$Sol[, l])
          tmp0 <- tmp0[round(seq(1, length(tmp0), by = thin))]
          #alternative: tmp0 <- tmp0[round(seq(1, length(tmp0), length = desiredlength))]
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
          tmp0 <- tmp0[round(seq(1, length(tmp0), by = thin))]
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
}
