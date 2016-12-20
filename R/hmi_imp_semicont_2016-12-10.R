#' The function for hierarchical imputation of semicontinous variables.
#'
#' The function is called by the wrapper.
#' @param y_imp_multi A Vector with the variable to impute.
#' @param X_imp_multi A data.frame with the fixed effects variables.
#' @param heap A numeric saying to which (single) values the data might be heaped.
#' @param M An integer defining the number of imputations that should be made.
#' @return A n x M matrix. Each column is one of M imputed y-variables.
imp_semicont<- function(y_imp_multi,
                             X_imp_multi,
                             heap = 0,
                             M = 10){


  tmp_data <- cbind(y_imp_multi, X_imp_multi)
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

  what_method <- ind.imp <- imp_binary(y_imp_multi = y_binary,
                                             X_imp_multi = X_imp_multi,
                                             M = M)
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
    y1_imp <-  imp_cont(y_imp_multi = y_imp_multi[what_method[, i] == 1],
                              X_imp_multi = X_imp_multi[what_method[, i] == 1, ],
                              M = 1)


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

