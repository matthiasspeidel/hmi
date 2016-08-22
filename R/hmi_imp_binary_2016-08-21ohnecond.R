### Funktion zur Anonymisierung/Imputation binaerer Variablen


################################################################################
### fehlende Werte in den Variablen werden zunaechst
### durch Single Imputation ersetzt.
################################################################################

################################################################################
## Argumente:
## y.variable.name: Name der zu anonymisierende Zielgroesse
## data.org.bin.imp: die zur Imputation verwendeten Daten (inkl der zu imputierenden Variable)
## cond: ??? Logische Variable die die Imputationsfunktion zwingt fuer Werte
## 0 oder 1 zu imputieren (???)
## mn: Mindesbesetzung der Kategorien von factor-Variablen. Faktor-Stufen mit
## einer absoluten Haeufigkeit < mn werden zu anderen Stufen zugewiesen
## max.se: Der erlaubte Standardfehler fuer die Schaetzung eines Regressionsparameters
## Wenn der Standardfehler zu dem Koeffizient einer Variable > max.se, wird diese
## Variable aus dem Imputationsmodell gestrichen

#' The function for non-hierarchical imputation of binary variables.
#'
#' This function is called by the wrapper.
#' @param y_variable_name A character specifying the variable to impute.
#' @param data_org_bin_imp A data.frame of the data to impute
#' (including the one with the missing values on the covariates.)
#' @param M An integer defining the number of imputations that should be made.
#' @param mn An iteger specifying the number of observation in each category of factor variables
#' QUESTION TO MY SELF: IS THIS THE RIGHT PLACE FOR THIS???
#' @param max_se A numeric specifying the maximal standard error of a covariate's parameter.
#' Covariates with a standard error > max.se are removed from the model for stability.
#' @return A n x M matrix. Each column is one of M imputed y-variables.
imp_binary <- function(y_variable_name, data_org_bin_imp,
                       M = 10, mn = 6, max_se = NULL){

  n <- nrow(data_org_bin_imp)
  #get y from the data
  y_org <- subset(data_org_bin_imp, select = y_variable_name)

  y_ret <- y_org

  #the missing indactor indicates, which values of y are missing.
  mis_indicator_y_org <- is.na(y_org)
  X_org <- subset(data_org_bin_imp, select = -which(names(data_org_bin_imp) == y_variable_name))

  #MS: schaue welche Variablen ein strukturelles/Filterbedingtes NA haben.
  #Bsp: Alter vom Kind = NA, wenn Anzahlkinder = 0.
  mis_indicator_X_org <- apply(X_org, 2, FUN = function(x) any(is.na(x)))

  #MS: entferne diese Variablen
  X_org <- subset(X_org, select = !mis_indicator_X_org)


  # Laden von Paketen

  if (requireNamespace("boot", quietly = TRUE)) { # for inv.logit
    inv.logit <- boot::inv.logit
  } else {
    inv.logit <- function(x) exp(x)/(1 + exp(x))
    print("We suggest to install the package 'boot'.")
  }

  requireNamespace("MASS", quietly = TRUE) # for mvrnorm

  ################################################################################
  ### Imputation
  ################################################################################

  y_obs <- y_org[!is.na(y_org), , drop = FALSE]
    #ALT: subset(y_org, bedingung_obs)
  #Hinweis: y.obs hat (fehlende Werte vorausgesetz) weniger Beobachtungen als y.cond

  #erste Modellmatrix f?r das Imputationsmodell getrennt nach
  #West/Ost und Betriebsgr??enquantil
  X_obs <- subset(X_org, !is.na(y_org)) # exkl. unbeobachteter Subjekte

  ## Kontrolle der Besetzung und ggf. Fusion von Faktorlevels
  ## Faktorvariablen in Regressormatrix
  cl <- which(lapply(X_obs, class) == "factor")


  ## Zugriff auf Faktoren
  #durchlaufe jeder Faktorvariable
  #!!!dauert relativ lang. evtl mit lapply etc. arbeiten!!!
  for(l3 in cl){

        #Haeufigkeiten der Faktorstufen-Auspraegungen bestimmen
        tbl0 <- table(X_obs[, l3])
        #Ermitteln der Anzahl der Faktorstufen, deren beobachtete Haeufigkeit
        #zu gering (kleiner als mn) ist
        d0 <- sum(tbl0 < mn)

        #die Levels mit zu wenigen Beobachtungen werden zusammengefasst
        levels(X_obs[, l3])[tbl0 < mn] = "my.other"
        #Note: es kann nat?rlich sein, dass "my.other" dann die gr??te Gruppe bildet...

        #wenn nun immer noch zu wenige Beobachtungen drin sind...,
        #Haeufigkeiten der Faktorstufen-Auspraegungen,
        #nach Erstellung von my.other bestimmen
        tbl1 <- table(X_obs[, l3])

        d1 <- sum(tbl1 < mn)
        if(d1 > 0){
          #... fuehre sie mit der kleinsten, gueltigen Auspr?gung zusammen.
          #make new temporary table without the my.other category
          tbl2 <- table(X_obs[, l3])[names(table(X_obs[, l3])) != "my.other"]

          #get the smallest (but valid) category
          smallest_cat <- which(tbl2 == min(tbl2))

          #rename the smallest categeory's level by the old levelname
          #combined with "my.other"  (STEP B)
          #and "my.other" get the same name
          new_name <- paste(levels(X_obs[, l3])[smallest_cat], "my.other", sep = "_")
          levels(X_obs[, l3])[levels(X_obs[, l3]) == "my.other"] <- new.name  #STEP A
          levels(X_obs[, l3])[smallest_cat] <- new_name   #STEP B
        }

        if(d0 > 0){
          print(paste(d0, "Faktorlevel(s) von", attributes(X_obs)$names[l3],
                      "mit Besetzung <", mn, "werden fusioniert."))
        }
  }

  #hier muesste man sich eigentlich auf die beobachteten Daten beschraenken,
  #denn angenommen ich habe folgenden Fall:
  #alle Personen, bei denen y beobachtet ist, haben x = "A",
  #die mit nicht beobachtetem y haben x = "B".
  #dann wuerde no.var auf den vollen Datensatz die Variable als nicht singulaer sehen.
  #Die Daten, die zum schaetzen des Modells genommen werden, sind aber singulaer!

  number_unique_values <- sapply(X_obs, function(x) length(unique(x)))
  #Funktion zur Erkennung singulaerer Variablen (also Variablen, bei denen
  #nur genau eine Merkmalsauspraegung beobachtet wurde,
  #und die somit keinerlei Erklaerungkraft fuer die Zielvariable haben kann)

  print(paste(sum(number_unique_values == 1),
              "Variable(s) were removed from the model because they were constant."))


  ## zweite Modellmatrix fuer das Imputationsmodell
  estimation_X0 <- subset(X_obs, select = which(number_unique_values > 1))
  prediction_X0 <- subset(X_org, select = which(number_unique_values > 1))
  # but now, add an Intercept for a better model fit.
  estimation_X0$Intercept <- 1
  prediction_X0$Intercept <- 1
  #rechne logistisches Modell bei dem y1
  #(y, die  wo == l1 & eval(cond) & mis.indicator == FALSE erf?llen)
  #durch X.reg (also alle X, die die selbe Bedingung erf?llen) erkl?rt wird.
  #anders formuliert: im ersten Durchlauf die responder aus dem Westen.
  #Im zweiten die Responder aus dem Osten.
  estimation_data0 <- data.frame(y_obs, estimation_X0)

  prediction_y <- sample_imp(y_org)
  prediction_data0 <- data.frame(prediction_y, prediction_X0)#data_org_bin_imp

  tmp_formula <- as.formula(paste(y_variable_name, "~ 0 + ."))

  #Schaetzung eines erstens Modells. Da Anpassungen wie Entfernen von Variablen mit
  #zu grossem Standardfehler vorgenommen werden, ist spaeter eine erneute
  #Modellschaetzung mit den angepassten Daten notwendig.
  estimation_model1 <- glm(tmp_formula, data = estimation_data0,
                           family = binomial(link = "logit"))

  prediction_model1 <- glm(tmp_formula, data = prediction_data0,
                          family = binomial(link = "logit"))

  estimation_X1 <- model.matrix(estimation_model1) # Note, that now an intercept is present
  #HINWEIS: MODELL WIRD NUR GESCHAETZT, UM DIE model.matrix ZU BEKOMMEN!
  # DIE PARAMETER DES MODELS SIND TOTAL EGAL UND NUTZLOS!
  prediction_X1 <- model.matrix(prediction_model1)

  #!!!damit man eine Model.matrix fuer die Individuen mit fehlendem y-Wert bekommt,
  #werden die Ypsilons nun einfach durch
  #Zufalsszahlen ersetzt.
  #zwischen tmp und X.imp1 kann sich die Anzahl der Variablen erhoehen,
  #weil die Faktor-Stufen als Dummy-Variablen in der der ModelMatrix
  #gespeichert werden

  ## Ueberpruefung ob NAs im Vektor geschaetzter Koeffizienten
  na_test <- is.na(estimation_model1$coeff)

  print(paste(sum(na_test),
              "Variable(s) were removed from the model because their coefficients were NA"))

  # dritte Modellmatrizen
  ## Variablen deren Koeffizient nicht gesch?tzt werden kann, werden entfernt
  ## (-1 wegen Intercept!!)
  prediction_X2 <- subset(prediction_X1, select = which(!na_test))
  estimation_X2 <- subset(estimation_X1, select = which(!na_test))

  estimation_data2 <- data.frame(y_obs, estimation_X2)

  prediction_data2 <- data.frame(prediction_y, prediction_X2)


  estimation_model3 <- glm(tmp_formula, data = estimation_data2,
                           family = binomial(link = "logit"))

  prediction_model3 <- glm(tmp_formula, data = prediction_data2,
                           family = binomial(link = "logit"))

  estimation_X3 <- model.matrix(estimation_model3)
  prediction_X3 <- model.matrix(prediction_model3)


  ### Kontrolle der Standardfehler der Parameter-Schaetzungen
  coef_std <- sqrt(diag(vcov(estimation_model3)))
  ## Entfernung der Variablen, deren Koeffizientenschaetzungen zu gross sind.
	#"zu gross" soll hier bedeuten: der geschaetzte Standardfehler ist 20 mal
  # (Zahl von Joerg) groesser als der geschaetzte Paramter
	#oder wenn max.se spezifiziert wurde, dann wird dieser fixe Wert herangezogen

  if(is.null(max_se)){
        max_se <- abs(coef(estimation_model3) * 1)
  }

  print(paste(sum(coef_std > max_se),
              "Variable(s) were removed from the model",
              "because their standard errors were too high."))

  estimation_X4 <- subset(estimation_X3, select = which(coef_std <= max_se))
  prediction_X4 <- subset(prediction_X3, select = which(coef_std <= max_se))

  ## Es sollen mind. doppelt so viele Beobachtungen vorhanden sein,
  #wie man Modellkoeffizienten sch?tzen m?chte
  ## ansonsten soll das Modell gar nicht erst gesch?tzt werden.
  if(!(nrow(y_obs) * 2 > ncol(estimation_X4))){
        stop("Too few observations in comparison to the number of covariates.
             Consider removing unpredictive variables.")
  }

  estimation_data4 <- data.frame(y_obs, estimation_X4)

  prediction_data4 <- data.frame(prediction_y, prediction_X4)

  estimation_model5 <- glm(tmp_formula, data = estimation_data4,
                           family = binomial(link = "logit"))

  prediction_model5 <- glm(tmp_formula, data = prediction_data4,
                           family = binomial(link = "logit"))

  prediction_X5 <- model.matrix(prediction_model5)


  ## geschaetze Regressionsparameter:
  reg_coef <- coefficients(estimation_model5)
  coef_var <- vcov(estimation_model5)

  y_imp <- array(NA, dim = c(n, M))
  ###start imputation
  for (j in 1:M){

    ###  Ziehung der Parameter aus ihrer Posterior Distribution

    ## Ziehung simulierter Regressionsparameter
    newbeta <- mvrnorm(1, reg_coef, coef_var)

    ## Berechnung von Wahrscheinlichkeiten f?r Y = 0 bzw. Y = 1
    #anhand der logistischen Gleichung
    #Prob(y = 1) = P
    P <- inv.logit(prediction_X5 %*% newbeta)

    u <- runif(n = n) #ziehe gleichverteilte Zufallsvariablen
    #Wenn nun die Berechnete Wahrscheinlichkeit groesser ist,
    #als die gleichverteilte Zufallsvariable, wird Y = 1, ansonsten 0.

    y_tmp <- ifelse(P >= u, 1, 0)

    # If our y.org are real factor levels like "male" and "female"
    # transform the 0s and 1s to those factors:
    if(is.factor(y_org[, 1])){
        y_tmp <- levels(y_org[,1])[y_tmp + 1]#!!!checken, ob das hier richtig ist, oder
        #ob es nicht 2-y_tmp heissen muesste!!!
    }

    y_imp[, j] <- ifelse(is.na(y_org), y_tmp, y_org[,1])

  }
  return(y_imp)
}

