#' The function for hierarchical imputation of semicontinous variables.
#'
#' The function is called by the wrapper.
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
#' @param allowed_max_value A single numeric Value which shall not be exceeded
#' when values are imputed (e.g. the age of a person can be limited to 125).
#' @param allowed_max_variable A character naming a variable V.
#' For each Y_i the value of V_i shall not exceeded
#' (e.g. the net income shall not exceed the gross income).
#' Note that a new imputed value has to satisfy both conditions of \code{allowed_max_value}
#' and \code{allowed_max_variable} at the same time.
#' @param allowed_min_value Analog to \code{allowed_max_value}.
#' @param allowed_min_variable Analog to \code{allowed_max_variable}.
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
                             burnin = 1000,
                             allowed_max_value = Inf,
                             allowed_max_variable = NULL,
                             allowed_min_value = -Inf,
                             allowed_min_variable = NULL){


  tmp_data <- cbind(y_imp_multi, X_imp_multi, Z_imp_multi, clID)
  n <- nrow(tmp_data)

  #the missing indactor indicates, which values of y are missing.
  mis_indicator <- is.na(y_imp_multi)
  #get the defaults values for heap
  if(is.null(heap)) heap = 0


  #these steps are neccesary, because in wrapper there is a value given for heap and max.se
  #but those values could be NULL


  #mache aus y.org eine binaere Variable, wobei 0 fuer heap; und 1 fuer ungleich heap steht
  y_binary <- y_imp_multi

  #Die Beobachtungen, die gleich heap und nicht fehlend sind...
  condition0 <- (y_imp_multi == heap) & !is.na(y_imp_multi)
  #...werden 0 gesetzt.
  y_binary[condition0] <- 0

  #Die Beobachtungen, die ungleich heap und nicht fehlend sind...
  condition1 <- (y_imp_multi != heap) & !is.na(y_imp_multi)
  #...werden 1 gesetzt.
  y_binary[condition1] <- 1

  #Use the imputation function of the binary variable on the indicator
  #to set ind.imp to 0 or 1
  #Zuerst ermitteln, ob 0 oder ein stetiger Wert imputiert werden soll.

  what_method <- ind.imp <- imp_binary_multi(y_imp_multi = y_binary,
                                             X_imp_multi = X_imp_multi,
                                             Z_imp_multi = Z_imp_multi,
                                             clID = clID,
                                             M = M,
                                             nitt = nitt,
                                             thin = thin,
                                             burnin = burnin,
                                             allowed_max_value = allowed_max_value,
                                             allowed_max_variable = allowed_max_variable,
                                             allowed_min_value = allowed_min_value,
                                             allowed_min_variable = allowed_min_variable)
  y_imp <- array(NA, dim = c(n, M))
  for(i in 1:M){


    #Es werden nun die Verwendet, die in der indikator-Imputation eine 1 bekommen haben.
    #(das sind die jenigen, deren wahrer Wert vermutlich nicht gleich heap ist)
    #In anderen Worten: alle die nach der Imputation 0 haben
    #(also in der original-Variable heap), werden aus der cont-imp ausgeschlossen.

    y_tmp <- what_method[, i]

  #decide whether I want to make the functions to give the full data sor the X, Y, Z, seperately

    # use the imputation function of the continuous variable to generate y1.imp
    #data.step1
    y1_imp <-  imp_cont_multi(y_imp_multi = y_imp_multi[what_method[, i] == 1],
                              X_imp_multi = X_imp_multi[what_method[, i] == 1, ],
                              Z_imp_multi = Z_imp_multi[what_method[, i] == 1, ],
                              clID = clID[what_method[, i] == 1],
                              M = 1,
                              nitt = nitt,
                              thin = thin,
                              burnin = burnin,
                              allowed_max_value = allowed_max_value,
                              allowed_max_variable = allowed_max_variable,
                              allowed_min_value = allowed_min_value,
                              allowed_min_variable = allowed_min_variable)


    # set the final value of y
    #Dabei nehme ich ind.imp zur Grundlage. Die imputierten (oder originalen) Nuller
    #koennen beibehalten werden. Die Einser hingegen muessen noch durch die Imputierten
    #(oder originalen) Werte ungleich 0 ersetzt werden.
    #Bei heap != 0, muessen natuerlich die NUller durch heap ersetzt werden


    y_tmp[what_method[, i] == 1] <- y1_imp

    y_imp[ , i] <- y_tmp

  }

  return(y_imp)
}

