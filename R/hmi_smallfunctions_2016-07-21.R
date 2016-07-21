
##############################
###Gibbs sampler functions####
##############################


####Gibbs sampler for general random coefficients models

#' helper function to be called in \code{update.alpha.coef}.
#'
#' It (basically) calculates (Z'Z + sigma.e * sigma.b^-1)^-1 * z' resid.
#' @param x A data.frame which is basically cbind(Z_obs, y_obs - X_obs * beta.new).
#' @param sigma.y.new The numeric residual variance.
#' @param sigma.alpha.new.inv The inverse random effects covariance matrix.
#' @return It returns (for this given cluster) a (cluster specific) random intercept and random slope.
calc <- function(x, sigma.y.new, sigma.alpha.new.inv){
  solve(t(x[, 1:(ncol(x) - 1)]) %*% as.matrix(x[, 1:(ncol(x) - 1)]) +
          sigma.y.new^2 * sigma.alpha.new.inv) %*% t(x[, 1:(ncol(x) - 1)]) %*% x[, ncol(x)]
}


#' update.alpha.coef
#'
#' The function updates the cluster specific random effects
#' by drawing from a multivariate normal distribution
#' (with a unique distribution for every cluster).
#' @param Z_obs The random effects data matrix of those observations with an observed value of y.
#' @param y_obs The target variable of those observations with an observed value of y.
#' @param X_obs The fixed effects data matrix of those observations with an observed value of y.
#' @param beta.new The vector of fixed effects parameters.
#' @param clID_obs The cluster ID vector of those observations with an observed value of y.
#' @param sigma.y.new The numeric residual variance.
#' @param sigma.alpha.new.inv The inverse random effects covariance matrix.
#' @return It returns cluster specific random effects.
update.alpha.coef <- function(Z_obs, y_obs, X_obs, beta.new, clID_obs, sigma.y.new, sigma.alpha.new.inv){

  # -- calculate cluster specific random effects
  tmp <- cbind(Z_obs, y_obs - X_obs %*% beta.new)
  alpha.hat <- by(tmp, clID_obs, calc,
                  sigma.y.new = sigma.y.new, sigma.alpha.new.inv = sigma.alpha.new.inv)

  rm(tmp)

  # -- calculate cluster specific random effects covariance matrices.
  var.alpha.hat <- by(Z_obs, clID_obs,
                      function(x) sigma.y.new^2 *
                        solve(t(x) %*% as.matrix(x) + sigma.y.new^2 * sigma.alpha.new.inv))

  # -- draw for each cluster new random effects.
  alpha.new <- t(mapply(MASS::mvrnorm, mu = alpha.hat, Sigma = var.alpha.hat, n = 1))
  return(alpha.new)
}



#' update.beta.coef
#'
#' The functions updates the fixed effects (imputation) parameters.
#' By drawing them from their normal distribution.
#' @param Z_obs A data.frame with the random effects variables.
#' @param alpha.new A matrix (with the dimension: number of clusters times number of random effect variables)
#' of cluster specific random effects.
#' @param b.hat.part1 The (X'X)^{-1} X' matrix.
#' @param clID_obs The cluster ID vector of those observations with an observed value of y.
#' @param y_obs The target variable of those observations with an observed value of y.
#' @param sigma.y.new The residual variance (positive numeric).
#' @param xtx The (X'X)^{-1} matrix.
#' @return It returns a vector of new fixed effects parameters.
update.beta.coef <- function(Z_obs, alpha.new, b.hat.part1, clID_obs, y_obs, sigma.y.new, xtx){
  Z.alpha <- rowSums(Z_obs * alpha.new[clID_obs,])
  beta.hat <- b.hat.part1 %*% (y_obs - Z.alpha)
  Sigma.hat <- sigma.y.new^2 * xtx
  newbeta <- MASS::mvrnorm(n = 1, mu = beta.hat, Sigma = Sigma.hat)
  return(newbeta)
}

#' update.sigma.y.coef
#'
#' The function updates the residual variance parameter by drawing from a chisquared distribution.
#' @param y_obs The target variable of those observations with an observed value of y.
#' @param X_obs The fixed effects data matrix of those observations with an observed value of y.
#' @param beta.new The vector of fixed effects parameters.
#' @param Z_obs The random effects data matrix of those observations with an observed value of y.
#' @param alpha.new The matrix of cluster specific random effects.
#' @param clID_obs The cluster ID vector of those observations with an observed value of y.
#' @param n.obs The number of individuals with an observed y.
#' @return The numeric residual variance.
update.sigma.y.coef <- function(y_obs, X_obs, beta.new, Z_obs, alpha.new, clID_obs, n.obs){
  resid <- y_obs - X_obs %*% beta.new - rowSums(Z_obs * alpha.new[clID_obs, ])
  sigma.y.new <- sqrt(sum(resid^2)/rchisq(1, n.obs - 1))
  return(sigma.y.new)
}

#' update.sigma.alpha.coef
#'
#' It updates the random effects covarariance matrix by drawing from a Wishart distribution.
#' @param hyper.df1 A vector of numeric hyper parameters
#' @param hyper.df2 A matrix of hyper parameters.
#' @param alpha.new A matrix of cluster specific random effects.
#' @return It returns a matrix with a new random effects covarariance matrix.
update.sigma.alpha.coef <- function(hyper.df1, hyper.df2, alpha.new){
  hyper.Sigma <- solve(t(alpha.new) %*% alpha.new + hyper.df2)
  new.sigma.alpha <- solve(rWishart(1, df = hyper.df1, Sigma = hyper.Sigma)[,,1])
  return(new.sigma.alpha)
}


#help function to do a sample imputation
#' Sample imputation.
#'
#' Function to sample values in a variable from other (observed) values in this variable.
#' So this imputation doesn't use further covariates.
#' @param variable A vector of size \code{n} with missing values.
#' @return A vector of size \code{n} without missing values.
#' @examples
#' sample_imp(c(1, NA, 3, NA, 5))
sample_imp <- function(variable){
  #!!!Es waere sinnvoll hier schon Nebenbedingungen zuzulassen!!!
  ret <- variable
  ret[is.na(ret)] <- sample(size = sum(is.na(variable)),
                             variable[!is.na(variable)], replace = TRUE)
  return(ret)
}


#TODO: Fallunterscheidungen um mit Vectoren, Matrizen oder data.frames umzugehen!!!
#Variable sollte ein Vektor oder factor sein

#' Get the type of variables.
#'
#' Function checks wether a variable is: ...
#' \itemize{
#'  \item continuous (just many numeric values),
#'  \item semicontinuous (many numeric values but more than 5% of the data are 0),
#'  \item a intercept (just the same value for all observations),
#'  \item binary (just two different values - like 0s and 1s or "m" and "f"),
#'  \item categorical (the variable is a factor or has more than 3 different values)
#'  \item orderd categorical (the categorical variable is ordered.)
#'}
#'
#' @param variable A variable from your dataset.
#' @return A character denoting the type of \code{variable}.
get_type <- function(variable){

  if(length(table(variable)) == 1){
    type <- "intercept"
    return(type)
  }

  if(length(table(variable)) == 2){
    type <- "binary"
    #!!!ggfs muessen binaere Variablen noch in 0 und 1 umgewandelt werden.
    #Bsp. wenn die Geschlechtervariable ist bis dahin "m" und "w" ist.
    return(type)
  }

  if(is.vector(variable) || is.factor(variable) || is.array(variable)){
    if(is.numeric(variable)){
      type <- "cont"

      #if the variable is numeric and more than 5% of the variable values are 0,
      #I count this variable as semi-continious.
      #???FROM WHICH RATE ON, ONE CAN SAY, THAT THE VARIABLE IS SEMI-CONT AND NOT
      #ONLY CONT? 1%, 10%, 50%???
      #!!!problematische Zahl frei w?hlbar!!!
      if(sum(variable == 0, na.rm = TRUE)/
         length(variable[!is.na(variable)]) > 0.05){
        type <- "semicont"
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
    #MS: Wenn es kein Vector oder kein Factor ist, gehe ich davon aus,
    #dass es eine Matrix oder ein data.frame ist.
    #print("R ist komisch!")
    #MS: d.h. ich schliesse aus, dass ein eine Liste mit Datensaetzen etc. ist.
    #MS: Fallunterscheidungen um mit Vectoren, Matrizen oder data.frames umzugehen:
    if(ncol(variable) == 1){
      ret <- get_type(variable[, 1])
      return(ret)
    }
    if(ncol(variable) > 1){
      if(is.matrix(variable)){
        variable <- as.data.frame(variable)
      }
      ret <- unlist(lapply(variable, get_type))
      #MS:?!? Fuer data.frames brauche ich lapply, fuer matrizen apply.
      return(ret)
    }
  }
}



#' Averages the results of the imputation function \code{wrapper}.
#'
#' This function applies the analysis the user is interested in, on all different imputed dataset.
#' Then the results are pooled by simply averaging the results. So the user has to make sure that
#' his analysis produces results with a meaningful average. And furthermore has to accept that no
#' variance is calculated for these parameters.
#' @param imputed_data A \code{mids} object from the \code{wrapper} imputation function.
#' @param analysis A user generated function that he is interested in. See examples.
#' @return A vector with all averaged results.
#' @examples
#' my.formula <- Reaction ~ Days + (1 + Days|Subject)
#' my_analysis <- function(complete_data){
#'  # In this list, you can write all the parameters you are interested in.
#'  # Those will be averaged.
#'  # So make sure that averaging makes sensen and that you only put in single numeric values.
#'  parameters_of_interest <- list()
#'
#'  # ---- write in the following lines, what you are interetest in to do with your complete_data
#'  my_model <- lmer(my.formula, data = complete_data)
#'
#'  parameters_of_interest[[1]] <- VarCorr(my_model)[[1]][1, 1]
#'  parameters_of_interest[[2]] <- VarCorr(my_model)[[1]][1, 2]
#'  parameters_of_interest[[3]] <- VarCorr(my_model)[[1]][2, 2]
#'  names(parameters_of_interest) <- c("sigma0", "sigma01", "sigma1")
#'
#'  # ---- do not write below this line.
#'  return(parameters_of_interest)
#'}
#'
#' test <- sleepstudy
#' test[sample(1:nrow(test), size = 20), "Reaction"] <- NA
#' hmi_imp <- wrapper(data = test, model_formula = my.formula)
#' hmifit <- hmi_pool(data = hmi_imp, analysis = my_analysis)
hmi_pool <- function(data, analysis){

  if (!is.mids(data)){
    stop("The data must have class mids.")
  }

  results <-list()

  for (i in 1:data$m) {
    results[[i]] <- analysis(complete(data, i))
  }

  tmp <- simplify2array(results)
  mode(tmp) <- "numeric"
  return(rowMeans(tmp))
}
