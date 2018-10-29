#' The function for hierarchical imputation of semicontinuous variables.
#'
#' The function is called by the wrapper. We consider data to be "semicontinuous" when
#' more than 5\% of the (non categorical) observations.\cr
#' For example in surveys a certain portion of people, when asked for their income,
#' report "0", which clearly violates the assumption of income to be (log-) normally distributed.
#' @param y_imp A Vector with the variable to impute.
#' @param X_imp A data.frame with the fixed effects variables.
#' @param Z_imp A data.frame with the random effects variables.
#' @param spike A numeric value saying to which values Y might be spiked
#' @param clID A vector with the cluster ID.
#' @param nitt An integer defining number of MCMC iterations (see MCMCglmm).
#' @param burnin burnin A numeric value between 0 and 1 for the desired percentage of
#' Gibbs samples that shall be regarded as burnin.
#' @param thin An integer to set the thinning interval range. If thin = 1,
#' every iteration of the Gibbs-sampling chain will be kept. For highly autocorrelated
#' chains, that are only examined by few iterations (say less than 1000).
#' @param pvalue A numeric between 0 and 1 denoting the threshold of p-values a variable in the imputation
#' model should not exceed. If they do, they are excluded from the imputation model.
#' @param k An integer defining the allowed maximum of levels in a factor covariate.
#' @return A list with 1. 'y_ret' the n x 1 data.frame with the original and imputed values.
#' 2. 'Sol' the Gibbs-samples for the fixed effects parameters.
#' 3. 'VCV' the Gibbs-samples for variance parameters.
#' @export
imp_semicont_multi <- function(y_imp,
                             X_imp,
                             Z_imp,
                             clID,
                             spike = NULL,
                             nitt = 22000,
                             burnin = 2000,
                             thin = 20,
                             pvalue = 0.2,
                             k = Inf){
  #If no spike was given, the list_of_spikes_maker is use, which basically returns the mode.
  if(is.null(spike)){
    spike <-  list_of_spikes_maker(data.frame(y_imp))$y_imp
  }

  #if for some near-to-impossible reasons, spike is still NULL, it is set to 0.
  if(is.null(spike)){
    spike <- 0
  }

  tmp_data <- cbind(y_imp, X_imp, Z_imp, clID)
  n <- nrow(tmp_data)

  #these steps are neccesary, because in wrapper there is a value given for spike and max.se
  #but those values could be NULL

  y_binary <- factor(rep(NA, length(y_imp)), levels = c(0, 1))

  #The observations that are equal to the spike and are not NA...
  condition0 <- (y_imp == spike) & !is.na(y_imp)
  #...are set to 0.
  y_binary[condition0] <- 0

  #The observations that are unequal to the spike and are not NA...
  condition1 <- (y_imp != spike) & !is.na(y_imp)
  #...are set to 1.
  y_binary[condition1] <- 1

  #Use the imputation function of the binary variable on the indicator
  #to set what_method to 0 or 1

  tmp1 <- imp_binary_multi(y_imp = y_binary,
                           X_imp = X_imp,
                           Z_imp = Z_imp,
                           clID = clID,
                           nitt = nitt,
                           burnin = burnin,
                           thin = thin,
                           pvalue = pvalue,
                           k = k)

  what_method <- tmp1$y_ret
  # use the imputation function of the continuous variable to generate y1.imp

  tmp2 <-  imp_cont_multi(y_imp = y_imp[what_method == 1],
                          X_imp = X_imp[what_method == 1, ,drop = FALSE],
                          Z_imp = Z_imp[what_method == 1, ,drop = FALSE],
                          clID = clID[what_method == 1],
                          nitt = nitt,
                          burnin = burnin,
                          thin = thin,
                          pvalue = pvalue,
                          k = k)

  y1_imp <- tmp2$y_ret

  # set the final value of y:
  # the observations with method 1 (continuous (non hepead) observation)
  # get the continuously imputed values
  y_tmp <- array(NA, dim = length(y_imp))
  y_tmp[what_method == 1] <- y1_imp[, 1]

  # the observations with method 0 (spiked observation)
  # get the spike
  y_tmp[what_method == 0] <- spike

  y_ret <- data.frame(y_ret = y_tmp)

  # --------- returning the imputed data --------------
  ret <- list(y_ret = y_ret, Sol = tmp2$xdraws, VCV = tmp2$variancedraws)
  return(ret)
}

