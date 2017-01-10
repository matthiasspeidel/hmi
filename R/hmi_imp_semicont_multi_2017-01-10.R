#' The function for hierarchical imputation of semicontinous variables.
#'
#' The function is called by the wrapper. We consider data to be "semicontinuous" when
#' more than 5\% of the (non categorical) observations.\cr
#' For example in surveys a certain portion of people, when asked for their income,
#' report "0", which clearly violates the assumption of income to be (log-) normally distributed.
#' @param y_imp_multi A Vector with the variable to impute.
#' @param X_imp_multi A data.frame with the fixed effects variables.
#' @param Z_imp_multi A data.frame with the random effects variables.
#' @param model_formula A \code{\link[stats]{formula}} used for the analysis model.
#' @param heap A numeric saying to which (single) values the data might be heaped.
#' @param clID A vector with the cluster ID.
#' @param M An integer defining the number of imputations that should be made.
#' @param nitt An integer defining number of MCMC iterations (see MCMCglmm).
#' @param thin An integer defining the thinning interval (see MCMCglmm).
#' @param burnin An integer defining the percentage of draws from the gibbs sampler
#' that should be discarded as burn in (see MCMCglmm).
#' @return A n x M matrix. Each column is one of M imputed y-variables.
imp_semicont_multi <- function(y_imp_multi,
                             X_imp_multi,
                             Z_imp_multi,
                             clID,
                             model_formula,
                             heap = 0,
                             M = 10,
                             nitt = 3000,
                             thin = 10,
                             burnin = 1000){


  tmp_data <- cbind(y_imp_multi, X_imp_multi, Z_imp_multi, clID)
  n <- nrow(tmp_data)

  #the missing indactor indicates, which values of y are missing.
  mis_indicator <- is.na(y_imp_multi)
  #get the defaults values for heap
  if(is.null(heap)) heap = 0


  #these steps are neccesary, because in wrapper there is a value given for heap and max.se
  #but those values could be NULL

  y_binary <- y_imp_multi

  #The observations that are equal to the heaping value and are not NA...
  condition0 <- (y_imp_multi == heap) & !is.na(y_imp_multi)
  #...are set to 0.
  y_binary[condition0] <- 0

  #The observations that are unequal to the heaping value and are not NA...
  condition1 <- (y_imp_multi != heap) & !is.na(y_imp_multi)
  #...are set to 1.
  y_binary[condition1] <- 1

  #Use the imputation function of the binary variable on the indicator
  #to set what_method to 0 or 1

  what_method <- imp_binary_multi(y_imp_multi = y_binary,
                                             X_imp_multi = X_imp_multi,
                                             Z_imp_multi = Z_imp_multi,
                                             clID = clID,
                                             M = M,
                                             nitt = nitt,
                                             thin = thin,
                                             burnin = burnin)
  y_imp <- array(NA, dim = c(n, M))
  for(i in 1:M){

    y_tmp <- what_method[, i]


    # use the imputation function of the continuous variable to generate y1.imp

    y1_imp <-  imp_cont_multi(y_imp_multi = y_imp_multi[what_method[, i] == 1],
                              X_imp_multi = X_imp_multi[what_method[, i] == 1, ,drop = FALSE],
                              Z_imp_multi = Z_imp_multi[what_method[, i] == 1, ,drop = FALSE],
                              clID = clID[what_method[, i] == 1],
                              M = 1,
                              nitt = nitt,
                              thin = thin,
                              burnin = burnin)


    # set the final value of y:
    # the observations with method 1 (continuous (non hepead) observation)
    # get the continuously imputed values
    y_tmp[what_method[, i] == 1] <- y1_imp

    y_imp[ , i] <- y_tmp

  }

  return(y_imp)
}

