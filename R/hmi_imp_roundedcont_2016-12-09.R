#Zu kl?ren: soll die Funktion mit Intervallangaben umgehen k?nnen?
#Dann br?uchte man entweder zwei variablen oder eine in der das intervallangaben und pr?zise
#d.h. intervalll?nge = 0 gemeinsam gespeichert sind.
#y <- c("1300", "2342", "[1000, 2000]", "2967", "NA") etc.

#Steinbruch
#imp_roundedcont_multi <- function(y.variable.name, data.org, intervall.obs, Y_lower, Y_upper, exclude_in_tmp.data2,
#                                  cond,
#                                  mn,
#                                  allowed.max.value,
#                                  allowed.max.variable,
#                                  allowed.min.value,
#                                  allowed.min.variable,
#                                  max.se = NULL, impsyn = "imp", MLestimator.output.path = NULL){


#' The function to impute rounded continuous variables
#'
#' For example the income in surveys is often reported rounded by the respondents.
#' See Drechsler, Kiesl and Speidel (2015) for more details.
#' @param y_imp_multi A Vector with the variable to impute.
#' @param X_imp_multi A data.frame with the fixed effects variables.
#' @param Z_imp_multi A data.frame with the random effects variables.
#' @param clID A vector with the cluster ID.
#' @param intercept_varname A character denoting the name of the intercept variable.
#' @param M An integer defining the number of imputations that should be made.
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
#' @return Currently a vector with the completed variable that was to be imputed
#' (so it includes the unchanged originaly observed data
#' and the values that have been missing, but now are imputed) !!!MAKE IT FIT TO THE OTHER IMPUATION ROUTINES!!!
imp_roundedcont <- function(y_imp_multi, X_imp_multi,
                                  intercept_varname = NULL, M,
                                  allowed_max_value = Inf,
                                  allowed_max_variable = NULL,
                                  allowed_min_value = -Inf,
                                  allowed_min_variable = NULL){

  #??????????????? MAYBE EXCLUDE CLUSTERS WITH TOO MANY CATEGORIES??????

  #######################################################################
  #MS: BEGIN get starting imputation values by maximizing the likelihood#

  missind <- is.na(y_imp_multi)


  types <- array(dim = ncol(X_imp_multi))
  for(j in 1:length(types)) types[j] <- get_type(X_imp_multi[, j])
  need_stand <- types == "cont"
  categorical <- types == "categorical"

  #remove categories with more than 10 observations as the model in the current form
  #will cause later numerical probles
  too_many_levels <- colnames(X_imp_multi[, categorical])[
    apply(X_imp_multi[, categorical], 2, function(x) nlevels(factor(x))) > 10]

  X_imp_multi <- X_imp_multi[, !names(X_imp_multi) %in% too_many_levels, drop = FALSE]

  X_imp_multi_stand <- X_imp_multi
  X_imp_multi_stand[, need_stand] <- scale(X_imp_multi[, need_stand])


  #MS: Einkommensvariable befuellen, damit mit ihr richtig gearbeitet werden kann...
  blob <- sample_imp(y_imp_multi)

  #MS: ...wie zum Beispiel eine Designmatrix erstellen.
  lmstart <- lm(blob ~ 0 + . , data = X_imp_multi_stand)

  MM_0 <- model.matrix(lmstart)
  #???WHY IS THERE catA, catB, catC and catD??? One of theme should be the reference category!!!
  #because i don't include an intercept!

  max.se <- abs(coef(lmstart) * 3)
  coef.std = sqrt(diag(vcov(lmstart)))


  includes_unimportants <- any(coef.std > max.se)
  counter <- 0
  while(includes_unimportants & counter <= ncol(X_imp_multi_stand)){
    counter <- counter + 1

    MM_0 <- as.data.frame(MM_0[, coef.std <= max.se, drop = FALSE])
    lm_less_variables <- lm(blob ~ 0 + . , data = MM_0)
    #remove regression parameters which have a very high standard error
    max.se <- abs(coef(lm_less_variables) * 3)
    coef.std = sqrt(diag(vcov(lm_less_variables)))

    includes_unimportants <- any(coef.std > max.se)
  }

  MM_1 <- as.data.frame(MM_0)



  inc <- y_imp_multi
  n <- length(y_imp_multi)
  mean.inc <- mean(inc, na.rm = TRUE)
  sd.inc <- sd(inc, na.rm = TRUE)
  inc.std <- (inc - mean.inc)/sd.inc

  log.inc <- log(inc)
  mean.log.inc <- mean(log.inc, na.rm = TRUE)
  sd.log.inc <- sd(log.inc, na.rm = TRUE)
  log.inc.std <- (log.inc - mean.log.inc)/sd.log.inc

  log.inc.std.tmp <- sample_imp(log.inc.std)#MS: wird nur benoetigt um die designmatrix zu kriegen (?)
  #lmstart2 <- lm(log.inc.std.tmp ~ 0 + ., data = MM_1)

  ##MS: preparing the ml estimation
  ###define rounding intervals

  round_base <- c(1, 5, 10, 50, 100, 500, 1000)
  intervals <- round_base/2

  #check if which observation are rounded
  #MS: %% berechnet den Modulo. Calculate the rounding degree only for those with not an missing value in inc


  #MS: New approach with a p that will be of length equal to nrow(final.data2)
  p1 <- y_imp_multi %% 5    ==  0  # divisable by 5
  p1[is.na(p1)] <- FALSE

  p2 <- y_imp_multi %% 10   ==  0  # divisable by 10
  p2[is.na(p2)] <- FALSE

  p3 <- y_imp_multi %% 50   ==  0  # etc
  p3[is.na(p3)] <- FALSE

  p4 <- y_imp_multi %% 100  ==  0  #
  p4[is.na(p4)] <- FALSE

  p5 <- y_imp_multi %% 500  ==  0  #
  p5[is.na(p5)] <- FALSE

  p6 <- y_imp_multi %% 1000 ==  0  #
  p6[is.na(p6)] <- FALSE

  p <- factor(p1 + p2 + p3 + p4 + p5 + p6, levels = c("0", "1", "2", "3", "4", "5", "6"))
  #MS: p ist Vektor der fuer jede Beobachtung (indirekt) angibt durch welchen Faktor sie ohne Rest teilbar ist.
  #MS: Bspw. bedeutet ein Wert von 4, dass der Wert von HEK0600 durch den Faktor 100 teilbar ist
  #MS: (weil er durch jeweils durch 5, 10, 50 und 100, teilbar ist. Also 4 mal wurde "TRUE" aufsummiert, was 4 ergbit.)
  ###indicator which variables need to be imputed #MS: because they are rounded (and not because they are missing)
  rounded <- p != 0
  #MS: !!! Dass sich rounded die Daten ohne intervall-Beobachtungen zur Grundlage hat, macht das spaetere Arbeiten etwas schwierig
  #MS: Q: Warum ist length(p > 0) == 8195, aber length(which(p > 0)) == 6980?
  #MS: A: Weil in which(p > 0) ja nur die Indizes der Elemente angegeben sind, die das Kriterium p > 0 erfuellen.
  #MS: vergleiche 1:10 >= 7 und which(1:10 >= 7)


  #####maximum likelihood estimation using starting values
  ####estimation of the parameters

  # estimation of the starting values for eta and the thresholds on the x-axis:
  # ordered probit maximum possible rounding on the rounded in income data

  probitstart <- MASS::polr(as.ordered(p[!missind]) ~ inc.std[!missind],
                      contrasts = NULL, Hess = TRUE, model = TRUE,
                      method = "probit")

  gammastart <- as.vector(probitstart$coefficients) # the fix effect(s)
  kstart <- as.vector(probitstart$zeta) # the tresholds (in the summary labeled "Intercepts")

  lmstart2 <- lm(log.inc.std[!missind] ~ 0 + ., data = MM_1[!missind, , drop = FALSE])
  betastart2 <- as.vector(lmstart2$coef)
  sigmastart2 <- summary(lmstart2)$sigma


  #####maximum likelihood estimation using the starting values
  #MS: Die ML-Schaetzer sind dann die Imputationsparameter fuer die improved imputation.

  function_generator <- function(para, X, y_in_negloglik, myp, mean.log.inc, sd.log.inc){
    ret <- function(para){
      ret_tmp <- negloglik2(para = para, X = X, y_in_negloglik = y_in_negloglik, myp = myp,
                              mean.log.inc = mean.log.inc, sd.log.inc = sd.log.inc)
      return(ret_tmp)
    }
    return(ret)
  }

  ###exclude obs below (above) the 0.5% (99.5%) income quantile before maximizing
  ###the likelihood. Reason: Some extrem outliers cause problems during the
  ###maximization
  #MS: bei der Imputation sollen sie spaeter vorhanden sein, weshalb hier ein Datensatz 'data.for.imp1' "abgespalten" werden koennte.
  #MS: Alternativ kann man diese Beobachtungen einfach bei der Maximierung (der Likelihood) ausschliessen,
  #MS: so dass nicht jeder weitere Schritt der Designmatrix-Erzeugung separat durchgefuehrt werden muss.
  #MS: Alle weiteren Schritte (vorallem die Desingmatrix-Erzeugung) werden deshalb fuer jeden Datensatz separat durchgefuehrt.

  quants <- quantile(y_imp_multi, c(0.005, 0.995), na.rm = TRUE)
  outliers <- which(y_imp_multi < quants[1] | y_imp_multi > quants[2])

  #MS: inc.intervall brauchen wir vermutlich hier nicht mehr.
  negloglik2_generated <- function_generator(para = para,
                                             X = MM_1[-outliers, , drop = FALSE],
                                             y_in_negloglik = y_imp_multi[-outliers],
                                             myp = as.numeric(as.character(p[-outliers])),
                                             mean.log.inc = mean.log.inc,
                                             sd.log.inc = sd.log.inc)
  #MS:!!!STARTWERTE NUR ZUM CHECKEN, OB TATSAECHLICH MAXIMUM GEFUNDEN, MIT 0.7 MULTIPLIZIERT!!!
  a <- Sys.time()
  print(paste("Starting ML-Maximization to derive imputation parameters", Sys.time()))
  m2 <- optim(par = c(kstart, betastart2, gammastart, sigmastart2), negloglik2_generated, method = "BFGS",
              control = list(maxit = 10000), hessian = TRUE)
  #MS: da optim() standardmaessig minimiert muessen wir die Likelihood-Funktion "umdrehen"
  #MS: also mit -1 multiplizieren.
  #MS und diese entstandene negative log-likelihood minimierem um das Maximum der Likelihood zu bekommen.
#if(!is.null(MLestimator.output.path)) save(m2, file = MLestimator.output.path)
  b <- Sys.time()
  print(paste("Time to find maximum likelihood estimates under improved model:",
              format(difftime(b,  a, units = "mins"))))
  #MS: dauert bei mir ca. 58 - 82 Minuten.
  par_ml2 <- m2$par
  hess <- m2$hessian


  #MS: Links zur nearest covariance matrix:
  #MS: http://quant.stackexchange.com/questions/2074/what-is-the-best-way-to-fix-a-covariance-matrix-that-is-not-positive-semi-defi
  #MS: nearPD(hess)$mat
  isSymmetric(Sigma_ml2)
  Sigma_ml2 <- diag(diag(solve(nearPD(hess)$mat), ncol = ncol(hess)))

  #det(hess)
  #M <- matrix(c(4, -2, -2, 3), ncol = 2)
  #solve(M)
  #M <- matrix(1:9, 3)
  #diag(M[seq(1, by = nrow(M) + 1, length = nrow(M))])

  #solve(diag(diag(M), ncol = ncol(M)))/det(M)

 # solve(diag(x = , ncol = 2))
  #solve(nearPD(hess, corr = TRUE)$mat) #!!!SOMETIMES NOT INVERTIBLE!!!


  #MS: END get starting imputation values by maximizing the likelihood#
  #####################################################################

  ###set starting values equal to the observed income
###rounded income will be replaced by imputations later
	inc.imp <- inc
	inc.std.imp <- inc.std
 	log.inc.std.imp <- log.inc.std

 	y.imp <- array(NA, dim = c(n, M))

 	for(j in 1:M){
 	  ####draw new parameters #MS: because it is a Bayesian imputation
    check <- TRUE
    #  counter <- 0
    while(check){
      pars <- mvtnorm::rmvnorm(1, mean = par_ml2, sigma = Sigma_ml2)
      #first eq on page 63 in Drechsler, Kiesl, Speidel (2015)
      #Can we work with a diagonal matrix as well, or is this too far from the posterior?

      ####test if drawn parameters for the thresholds are in increasing order
      ####and if the standard deviation of the residuals is<0
      ####if yes, draw again
      #MS: pars entspricht c(kstart, betastart2, gammastart, sigmastart2) (?)
      test <- c(pars[2:6] - pars[1:5], pars[length(pars)])
      #MS: pars[2:6] - pars[1:5] bewirkt, dass k1 - k0 und k2 - k1, ..., k5 - k4 gerechnet wird
      #MS: und das soll nicht kleiner 0 sein, weil k1-k0 < 0 aequivalent zu k1 < k0 und das geht ja nicht.
      #MS: Das letzte Element ist sigmastart2 was logischerweise positiv sein muss.
      check <- any(test < 0) #MS: Es muessen alle Werte in test >= 0 sein, sonst terminiert die Schleife nicht
      #    counter <- counter +1
      #    print(counter)
    }


    beta_hat <- as.matrix(pars[7:(length(pars) - 2)], ncol = 1)
    #MS: betastart2 hatte noch 33 Elemente, beta.hat hat nur noch 26 ?!?
    #MS: es stimmt aber, dass die ersten 6 Elemente zu kstart gehoeren und die letzten 2 sind eta und sigma.
    gamma1_hat <- pars[length(pars) - 1]
    sigma_hat <- pars[length(pars)]
    mu_g <- gamma1_hat * as.matrix(MM_1) %*% beta_hat
    mu_y <- as.matrix(MM_1) %*% beta_hat
    mymean <- cbind(mu_g, mu_y)

    #The covariance matrix from equation (3)
    Sigma <- matrix(c(1 + gamma1_hat^2 * sigma_hat^2,
                      gamma1_hat * sigma_hat^2, gamma1_hat * sigma_hat^2,
                      sigma_hat^2), nrow = 2)

    ###################################################################################
    #MS: BEGIN IMPUTATION ONLY FOR THOSE OBSERVATION WITH NO INCOME INFORMATION AT ALL#
    #MS: HERE NO UNROUNDING OR IMPUTING INTERVALL-DATA TAKES PLACE#####################

    #MS: Es werden zwar zu viele Imputationswerte gezogen, da viele von den Beobachtungen die HEK0600.miss
    #MS: haben auch Intervall-Angaben gemacht haben. Und diese Beobachtungen werden danach nochmals imputiert.

    mytry <- rnorm(n = sum(missind),
                   mean = as.matrix(MM_1[missind, , drop = FALSE]) %*% beta_hat, sd = sigma_hat)

    #MS: Vorschlags-Imputationswert:
    imp_temp <- exp(mytry * sd.log.inc + mean.log.inc)

    #MS: Einen Check, ob die gezogenen Imputierten Werte auch in dem tatsaechlich beobachteten Intervall liegen
    #MS: braucht man nicht, denn schliesslich zieht man ja bereits aus der trunktierten Verteilung.

    inc.imp[missind] <- imp_temp
  #MS: !!!setzt aktuell vorraus, dass inc.imp kein data.frame sondern ein Vektor ist!!!


    #MS: END IMPUTATION ONLY FOR THOSE OBSERVATION THAT NO INCOME INFORMATION AT ALL#
    #MS: HERE NO UNROUNDING OR IMPUTING INTERVALL-DATA TAKES PLACE####################
    ##################################################################################


    #################################
    #MS: BEGIN UNROUNDING-IMPUTATION#
    ###define bounds for the rounding basis
    bounds_hat <- c(-Inf, pars[1:6], Inf)
    ###define interval bounds for maximum possible rounding intervals
    #MS: Needed in the following imputation lopp
    #MS: but only for those observations who answered this question
    #MS: I cannot write 'log.inc - log(intervalls)' etc. because we calculate log(inv - intervals)
    #MS: und das kann nur nach 'log(inc) - log(1 - intervals/inc) augeloest werden.

    y_lower <- (log(y_imp_multi - intervals[as.numeric(as.character(p)) + 1]) -
                  mean.log.inc)/sd.log.inc

    y_upper <- (log(y_imp_multi + intervals[as.numeric(as.character(p)) + 1]) -
                  mean.log.inc)/sd.log.inc

    g_lower <- bounds_hat[as.numeric(as.character(p)) + 1]
    g_upper <- bounds_hat[as.numeric(as.character(p)) + 2]


    ###loop over all observations that need to be unrounded
    for (i in which(rounded)){
      #	if(k %% 1000 == 0){
      #			print(k)
      #		}

      test <- TRUE
      while(test){

        ###draw from truncated multivariate normal
        ###drawn y must be between y_lower and y_upper
        ###drawn g must be smaller than g_upper (g > g_upper is not consistent with
        ###rounding observed in the data)
        mytry <- tmvtnorm::rtmvnorm(1,
                          mean = mymean[i, ],
                          sigma = Sigma,
                          lower = c(-Inf, y_lower[i]),
                          upper = c(g_upper[i], y_upper[i]),
                          algorithm = "gibbs", burn.in.samples = 1000)

        ###draws with rejection sampling. If intervall restrictions are not fulfilled after 1 mio
        ###draws, NA is given back. In this case, the rounded data will stay unmodified
        if(is.na(mytry[1])){
          print(paste("corrected imputation not possible for record:", i))
          print(paste("observed income for that record:", y_imp_multi[i]))
          ##generate a rounding indicater that is always consistent with the observed data
          g_temp <- bounds_hat[2] - 1

          mytry <- c(g.temp, log.inc.std.imp[k])
        }
        ####get imputed rounding indicator
        round_int <- sum(mytry[1] > bounds_hat)

        ###get imputed income on original scale
        imp_temp <- exp(mytry[2] * sd.log.inc + mean.log.inc)

        ###test if imputed income after rounding is equal to observed rounded income given the
        ###imputed rounding indicator
        ###if not(test=TRUE), draw again
        #MS: test ist ein Indikator ob der imputierte Wert zur Annahme bezueglich des Rundungsprozess passt.
        test <- round(imp_temp/round_base[round_int]) * round_base[round_int] !=
          y_imp_multi[i]
      }

      inc.imp[i] <- imp_temp
    }
    y.imp[, j] <- inc.imp
  }
  #MS: END UNROUNDING-IMPUTATION#
  ###############################
  #MS: replace original income by imputed income
  return(y.imp) #cbind(inc.imp, final.data2[, -which(names(final.data2) == y.variable.name)])
}


