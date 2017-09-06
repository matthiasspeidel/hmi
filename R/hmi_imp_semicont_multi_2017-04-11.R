#' The function for hierarchical imputation of semicontinuous variables.
#'
#' The function is called by the wrapper. We consider data to be "semicontinuous" when
#' more than 5\% of the (non categorical) observations.\cr
#' For example in surveys a certain portion of people, when asked for their income,
#' report "0", which clearly violates the assumption of income to be (log-) normally distributed.
#' @param y_imp A Vector with the variable to impute.
#' @param X_imp A data.frame with the fixed effects variables.
#' @param Z_imp A data.frame with the random effects variables.
#' @param heap A numeric saying to which (single) values the data might be heaped.
#' @param clID A vector with the cluster ID.
#' @param nitt An integer defining number of MCMC iterations (see MCMCglmm).
#' @param thin An integer defining the thinning interval (see MCMCglmm).
#' @param burnin An integer defining the percentage of draws from the gibbs sampler
#' that should be discarded as burn in (see MCMCglmm).
#' @return A n x 1 data.frame with the original and imputed values.
imp_semicont_multi <- function(y_imp,
                             X_imp,
                             Z_imp,
                             clID,
                             heap = 0,
                             nitt = 3000,
                             thin = 10,
                             burnin = 1000){


  tmp_data <- cbind(y_imp, X_imp, Z_imp, clID)
  n <- nrow(tmp_data)

  #the missing indactor indicates, which values of y are missing.
  mis_indicator <- is.na(y_imp)
  #get the defaults values for heap
  if(is.null(heap)) heap = 0


  #these steps are neccesary, because in wrapper there is a value given for heap and max.se
  #but those values could be NULL

  y_binary <- y_imp

  #The observations that are equal to the heaping value and are not NA...
  condition0 <- (y_imp == heap) & !is.na(y_imp)
  #...are set to 0.
  y_binary[condition0] <- 0

  #The observations that are unequal to the heaping value and are not NA...
  condition1 <- (y_imp != heap) & !is.na(y_imp)
  #...are set to 1.
  y_binary[condition1] <- 1

  #Use the imputation function of the binary variable on the indicator
  #to set what_method to 0 or 1

  what_method <- imp_binary_multi(y_imp = y_binary,
                                  X_imp = X_imp,
                                  Z_imp = Z_imp,
                                  clID = clID,
                                  nitt = nitt,
                                  thin = thin,
                                  burnin = burnin)

  # use the imputation function of the continuous variable to generate y1.imp

  y1_imp <-  imp_cont_multi(y_imp = y_imp[what_method == 1],
                              X_imp = X_imp[what_method == 1, ,drop = FALSE],
                              Z_imp = Z_imp[what_method == 1, ,drop = FALSE],
                              clID = clID[what_method == 1],
                              nitt = nitt,
                              thin = thin,
                              burnin = burnin)


  # set the final value of y:
  # the observations with method 1 (continuous (non hepead) observation)
  # get the continuously imputed values
  y_tmp <- data.frame(what_method)
  y_tmp[what_method == 1, 1] <- y1_imp

  y_ret <- data.frame(y_imp = y_tmp)

  return(y_ret)
}

