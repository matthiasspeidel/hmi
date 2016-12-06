#!!!DIE FUNKTION ERZEUGT BISHER NUR EINEN DATENSATZ MIT IMPUTIERTEN WERTEN!!!
#!!!NACH KONVERGENZ MUESSTE MAN NOCH M * 10 oder 20
# (ODER WAS AUCH IMMER - JE NACH DEM WIE STARK DIE AUTOKORRELATION IST)
# DATENSAETZE ERZEUGEN UND JEDEN 10ten ODER 20ten DATENSATZ AUSGEBEN.

#Q: Man schaut sich die Autokorrelation von was an?
#Die Autokorrelation des Intercepts in jeder Variable fuer sich
# (so wie auch die Konvergenz ermittelt wird) und
#nimmt die maximale Schrittweite?
#Oder berechne ich ein Lineares Modell 'fake_y [Tilde] alle imputierten Variablen' und schau mir da
#die Autokorrelation des Intercepts an?

# Grundidee:
# Imputiere eine Variable nur mit vollstaendigen anderen Variablen
# Implementierung:
# vor der Iterationsschleife eine for-Schleife einbauen, die eine
# Liste mit vollstaendigen Variablen erstellt und alle Variablen nach Missing Rate (MR) sortiert.
# Wenn die Liste leer ist (d.h. zu Beginn haben alle Variablen fehlende Werte),
# dann ziehe ein Sample fuer die Variable mit der kleinsten MR.
# Danach wird die Liste aktualisiert. Anschliessend wird die Variable mit der
# zweitkleinsten MR imputiert. usw.

## parameter_list: Liste, mit Elementen die heissen wie die Datensatzvariablen,
# und wiederum jeweils eine Liste sind mit Elementen: cond, mn, max_se,
# allowed.max, allowed_max_variable, allowed.min und allowed_min_variable.
# Die Mindest- und Hoechstwerte beziehen sich aber nur auf die zu imputierenden Werte!
# beobachtete Werte > allowed_max_value bzw. < allowed_min_value werden nicht veraendert!
# Die Erfuellung solcher Bedingungen ist wie spaeter beschrieben, die Aufgabe des Nutzers
# vor dem Start der Imputation.
# Details zu den Parametern der Elemente von parameter_list:

## mn: Mindesbesetzung der Kategorien von factor-Variablen. Faktor-Stufen mit
## einer absoluten Haeufigkeit < mn werden zu anderen Stufen zugewiesen
## max_se: Der erlaubte Standardfehler fuer die Schaetzung eines Regressionsparameters
## Wenn der Standardfehler zu dem Koeffizient einer Variable > max_se, wird diese
## Variable aus dem Imputationsmodell gestrichen.
## (DAS SOLL NICHT SIGNIFIKANTE VARIABLEN AUS DEM IMPUTATIONMODELL FERNHALTEN!?)
## allowed_max_value: Fuer die Variable sollen pauschal keine Werte > allowed.max
## imputiert werden (bspw. kann man das Alter einer Person auf 150 Jahre beschraenken).
## Wurde in einem Schritt fuer eine (oder mehrere) Beobachtung(en) ein Wert > allowed.max
## als Impuationsvorschlag gezogen,
## so wird erneut gezogen, bis der Vorschlagswert <= allowed.max.
## allowed_max_variable: der Name einer Variable, die fuer jedes Individuum einen
## individuellen Hoechstwert als Impuationswert erlaubt.
## Beispiel: das Nettoeinkommen soll imputiert werden. Als allowed_max_variable
## wuerde sich dann das Bruttoeinkommen anbieten.
## die Vorgaben aus allowed_max_value und allowed_max_variable muessen dabei gleichzeitig erfuellt sein;
## also Vorschlagswert <= pmin(allowed_max_variable, allowed_max_value)
## allowed_min_value: analog zu allowed_max_value
## allowed_min_variable: analog zu allowed_max_variable


# Vor dem Aufruf von wrapper() sollte der Anwender
# 1 vollstaendig fehlende Beobachtungen entfernen.
# 2 Den Datensatz auf logische Konsistenz ueberpruefen
# 2.1 Unlogische Angaben NA setzen oder ggfs. die entsprechende Beobachtung
#    komplett aus dem Datensatz entfernen
# 2.2 Fehlende Angaben, die sich aus anderen Angaben zweifelsfrei rekonstruieren lassen
#    vervollstaendigen. Bspw. kann bei angegebener Schwangerschaft, die
#    fehlende Angabe zum Geschlecht mit "weiblich" ausgebessert werden
#    (sofern man der Angabe ueber die Schwangerschaft vertraut)


# 2016-07-12: Ich lege fest, dass das package zunächst so schlank wie moeglich sein soll.
# Deshalb werde ich 1. nur erlauben, dass man einen vollstaendigen data.frame + model_formula
# der funktion uebergibt.
# Und nicht wie bisher alternativ einzelen x, y, oder z variablen.
# 2. Faellt die Wahlmoeglichkeit zwischen Imputation und Synthetisierung weg.

#' The main warpper function called by the user. LATER THE USER WILL ONLY USE THE WRAPPER.
#'
#' The user has to passes to the function his data.
#' Optionally he pass his analysis model formula so that hmi runs the imputation model
#' in line with his analysis model formula.
#' And of course he can specify some parameters for the imputation routine
#' (like the number of imputations and iterations and the burn in percentage.)
#' The standard usage should be that the user gives his complete dataset
#' and his analysis model. But he also could just give y, X and Z and the cluser ID.
#' @param data A matrix or (better) a data.frame with all variables appearing in \code{model_formula}.
#' @param model_formula A \code{\link[stats]{formula}} used for the analysis model.
#' @param parameter_list A LIST OF PARAMETERS ADD MORE DETAILS.
#' cond ???Vielleicht ein consitenzcheck??? bspw. muss man den zustand
#' "kriegen die Kinder Kindergeld? JA=1/NEIN=0?" nur imputieren,
#' Wenn die Bedingung "Anzahl Kinder >= 1" erfüllt ist.
#' Frage: Wollen wir solche Konsitenzen anbieten? Wenn ja wie?
#' Ich wuerde dann cond in "impute_if" etc. umbenennen.
#' Es könnte bspw. der Nutzer entweder fest vorgeben, welche Beobachtungen imputiert werden sollen
#' Oder er gibt flexibel nur einen charakter an wie "daten$anzahlkinder >= 1" und dieser wird
#' dann jedesmal neu ausgewertet.
#' D.h. AUCH ALS ANFORDERUNG AN WRAPPER ALLGEMEIN:
#' DAS PROGRAMM MUSS MIT DESIGNBEDINGTEN NAs UMGEHEN KÖNNEN!
#' D:H: AUCH DASS IN ANDEREN VARIABLEN (ALS DER ZU IMPUTIERENDEN VARIABLE) NAs vorkommen können.
#' aktuell (2016-08-19) erscheint mir das eher als ein nice-to-have als ein must-have.
#' deshalb werde ich das zunächst nicht uebernehmen/implementieren.
#' mn An iteger specifying the number of observation in each category of factor variables
#' max.se A numeric specifying the maximal standard error of a covariate's parameter.
#' Covariates with a standard error > max.se are removed from the model for stability.
#' @param M An integer defining the number of imputations that should be made.
#' @param nitt An integer defining number of MCMC iterations (see MCMCglmm).
#' @param thin An integer defining the thinning interval (see MCMCglmm).
#' @param burnin An integer defining the percentage of draws from the gibbs sampler
#' that should be discarded as burn in (see MCMCglmm).
#' @return \code{M} data.frames. Each one has the same format as the original \code{data}.
#' But the difference is that in each data.frame the NAs are replaced by imputed values.
#' !!!THE RESULT SHOULD BE IN A COMMON FORMAT EASY TO COMBINE
#' (damit meine ich, dass man ohne großen Aufwand auf jeden Datensatz die inhaltliche Analyse anwenden kann
#' und die Ergebnisse dann kombiniert werden. Am besten schaue ich mir an, wie mice etc. das macht.
#' mice gibt ein \code{mids} object aus). ES IST UNTERSCHEIDET SICH ABER VON DEN MIDS OBJECTEN AUS mice,
#' DA WIR AUF MANCHE ELEMENTE (method) VERZICHTEN
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
#' pool(with(data = hmi_imp, expr = lmer(Reaction ~ Days + (1 + Days | Subject))))

wrapper <- function(data, model_formula = NULL, parameter_list = NULL,
                    M = 10,
                    nitt = 3000,
                    thin = 100,
                    burnin = 1000,
                    form = vector("character", length = ncol(data)),
                    post = vector("character", length = ncol(data))){

  my_data <- data #wird spaeter ausgegeben. Entspricht von der Form dem
  #original Datensatz, nur, dass alle fehlenden Werte (mit Ausnahme der design
  #bedingt fehlenden) durch imputierte Werte ersetzt wurden.

  #des weitern gibt es noch einen tmp_data Datensatz

  if(is.matrix(my_data)){
    warning("We need your data in a data.frame (and not a matrix).
            So we changed it automatically, but it would be better if you do it.")
    my_data <- as.data.frame(my_data)
  }

  #"Beobachtungen" bei denen alle Angaben fehlen, sollten im Sinne des Anwenders
  #entfernt werden. Da er aber vielleicht einen Grund hat sie zu uebergeben, wird
  #nur eine Warnung ausgegeben.
  if(any(apply(my_data, 1, function(x) all(is.na(x))))){
      warning("Some observations have a missing value in every variable.
              This can increase your estimate's variance.
              We strongly suggest to remove those observations.")
  }


  # Set up default values for key variables like the target variable, the clusterID,
  # and the random covariates
  multilevel_target_varname <- NULL
  #extract key variables from formula
  if(!is.null(model_formula)){

    # if it was a character string, make it a formula
    if(class(model_formula) == "character"){
      warning("We need your model_formula in a formula (and not a character).
            So we changed it automatically, but it would be better if you do it.")

      model_formula <- formula(model_formula)
    }

    if(!all(all.vars(terms(model_formula)) %in% names(my_data))){
      warning("Some variable(s) in model_formula not found in dataset.
           Please check your model_formula and data.")
    }

    # ----------Target variable
    multilevel_target_varname <- as.character(model_formula)[2]

    # check if this multilevel_target_varname realy occurs in the dataest
    if(!multilevel_target_varname %in% names(my_data)){
      stop(paste(target_varname, "not found in dataset. Please check your model_formula and data."))
    }


    # -------- Cluster ID
    clID_varname <- as.character(lme4::findbars(model_formula)[[1]])[3]
    #check if there is a cluster variable specified:...
    if(is.na(clID_varname)){
      #... if it is not, we dont have a random effects model
      clID_varname <- NULL
      random_intercept_exists <- FALSE
      randomeffects_varname <- NULL

    }else{
      # check if this variable name realy occurs in the dataest
      if(!(clID_varname %in% names(my_data))){
        stop(paste(clID_varname, "not found in dataset. Please check your model_formula and data."))
      }


      # ----------- Z
      #library(formula.tools) check lhs(model_formula) or r1 <- rhs.vars(model_formula)
      # and then r1[grep("\\|", r1)]
      ###MS: get the random effects variables
      #MS: split formula before the |
      # TODO: FUNKTION MUSS NOCH MIT KLARKOMMEN, WENN KEIN MODEL ANGEGEBEN WURDE
      #ODER WENN ES KEIN LMER MODELL IST!!!
      nearly_randomeffects_varname <-
        strsplit(as.character(lme4::findbars(model_formula)), "\\|")[[1]][1]
#model_formula und my.model etc. aufräumen!!!
      #MS: take all variable.names (which are seperated by +)
      #MS: later it might be worth to expand it to interactions etc.
      nearly_randomeffects_varname <- strsplit(nearly_randomeffects_varname, "\\+")[[1]]

      #MS: remove spaces
      randomeffects_varname <- gsub(" ", "", nearly_randomeffects_varname)

      # -- check for random intercept
      # If not explicitely removed a random intercept is always present
      # cf. lmer(mpg ~ cyl + (drat|am), data = mtcars)
      random_intercept_exists <- TRUE
      # Only if a 0 is in the model_formula, no random intercept is estimated
      # cf. lmer(mpg ~ cyl + (0 + drat|am), data = mtcars)

      if("0" %in% randomeffects_varname){
        random_intercept_exists <- FALSE
        #exclude random effects called "0"
        randomeffects_varname <- randomeffects_varname[randomeffects_varname != "0"]
      }

      # To be fool proof, if the user puts a 1 in his model_formula,
      # a random intercept shall be calculated, even if the model could do something different
      # eg in the case lmer(mpg ~ cyl + (1 + 0 + drat|am), data = mtcars)
      if("1" %in% randomeffects_varname) random_intercept_exists <- TRUE


      if(all(c("0", "1") %in% randomeffects_varname)){
        stop("Your model_formula includes 0 for no random intercept and 1 for a random intercept at the same time.
           Please decide whether you want a random intercept or not.")
      }

    }

    # ---------X variables

    fixedeffects_varname <- all.vars(delete.response(terms(lme4::nobars(model_formula))))
    # --check if intercept is present
    fixed_intercept_exists <- attributes(delete.response(terms(lme4::nobars(model_formula))))$intercept == 1

    # --remove cluster ID from the X-variables (but it shouldn't be there)
    fixedeffects_varname <- fixedeffects_varname[fixedeffects_varname != clID_varname]

    if("0" %in% fixedeffects_varname){
      fixed_intercept_exists <- FALSE
      #exclude fixed effects called "0"
      fixedeffects_varname <- fixedeffects_varname[fixedeffects_varname != "0"]
    }

    ## Handling of a intercept variable.
    # MS: for the imputation functions it is important,
    # that constant variables are handled correctly
    constant_variables <- apply(my_data, 2, function(x) length(unique(x)) == 1)

    if(sum(constant_variables) >= 2){
      stop(paste("Your data has two or more constant variables: ", names(my_data)[constant_variables],
                 ". One needs to be dropped to avoid multicolinearity.", sep = ""))
    }

    # make shure, the data will have a "Intercept" variable, if the model_formula says to model an intercept.
    intercept_varname <- NULL
    if(fixed_intercept_exists | random_intercept_exists | sum(constant_variables) == 1){

      # give intercept_varname a default value.
      intercept_varname <- "Intercept"
      #MS:?!? If there is already an intercept, consider renaming it to "(Intercept)"
      # for reasons of consitency with further code.
      #MS: If there is no intercept yet in the data, generate one

      if(sum(constant_variables) == 1){
        intercept_varname <- names(my_data)[constant_variables]
        print(paste("We interprete", intercept_varname,
                    "as the intercept variable and set its value to 1."))

        if((any(fixedeffects_varname == "1") & names(my_data)[constant_variables] != "1")|
           (any(randomeffects_varname == "1") & names(my_data)[constant_variables] != "1")){
          warning(paste("Your model_formula includes '1' for an intercept.
                        The constant variable in your data is named ", names(my_data)[constant_variables], ".
                        Consider naming the variable also '1'.", sep = ""))
        }

      }

      #make shure that the random effects variable exist, and has the value 1
      my_data[, intercept_varname] <- 1

      if(fixed_intercept_exists){
        # MS: man muss sicher stellen, dass der alte intercept name nicht noch zusätzlich auftaucht
        # bspw war bei den random effects "1" ein Element von randomeffects_varname.
        # der Intercept war aber "Intercept", so dass die Schleife Stand 2016-07-17
        # "Intercept" "1" "X"  ausgab.
        fixedeffects_varname <- c(intercept_varname, fixedeffects_varname)

        # If the intercept_varname is not "1",
        # then every fixedeffects_varname that is "1" has to be removed
        if(intercept_varname != "1"){
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


      # Replace the "1" in randomeffects_varname by the name of the intercept variable.
      fixedeffects_varname[fixedeffects_varname == "1"] <- intercept_varname
      fixedeffects_varname <- unique(fixedeffects_varname)

      # Replace the "1" in randomeffects_varname by the name of the intercept variable.
      randomeffects_varname[randomeffects_varname == "1"] <- intercept_varname
      randomeffects_varname <- unique(randomeffects_varname)

    }


    #check if this variable name realy occurs in the dataset
    if(!all(randomeffects_varname %in% names(my_data))){
      warning(paste("We didn't find",
                   paste(randomeffects_varname[!randomeffects_varname %in% names(my_data)],
                         collapse = " and "),
                   "in your data. Please check your model_formula and data."))
    }


    #MS: now safe the variables used for modeling fixed effects
    X <- my_data[, fixedeffects_varname, drop = FALSE]

    #MS: now safe the variables used for modeling random effects
    Z <- my_data[, randomeffects_varname, drop = FALSE]


  }else{ # wenn model_formula NULL ist.
	  clID_varname <- NULL
  }


#!!!hier muesste man ueberlegen, ob man schon jetzt alle variablen mit einem design
  #bedingtem Fehlen faktorisiert. Das sind die Variablen, bei denen in der
  #parameter list cond nicht NULL ist.
  #Man muesste dann halt wieder aufpassen, wenn diese variable imputiert werden soll
  #denn dann sollen die stetigen informationen nicht laenger kategorisiert sein.
#neue Variablen typen "filtered.cont", "filtered.binary", "filtered.semicont"
  #und "filtered.categorical" koennten sinnvoll sein.
#

  #get missing rates
  missing_rates <- apply(my_data, 2, function(x) sum(is.na(x)))/nrow(my_data)
  #MS: it's not necessary to divide by nrow(my_data)

  #get variables with missings
  incomplete_variables <- names(my_data)[missing_rates > 0]

  #MS: Nonresponse-Matrix aufbauen: (wird später für die Convergenz Diagnostic gebraucht)
  nonresponse_matrix <- is.na(my_data)
  #if all missing rates are > 0...
  if(all(missing_rates > 0)){
    #...make a sample imputation for the variable with the lowest missing rate
    #If two or more variables have the same minimal MR,
    #sample wich variable should be imputed with sample.imp
    tmp_variable <- sample(size = 1,
                           names(my_data)[which.min(missing_rates)])

    #now do the sample.imp
    my_data[, tmp_variable] <- sample_imp(variable = my_data[, tmp_variable])
    #?!? ueberlegen, ob man hier schon sofern moeglich einen Plausibilitaetscheckt
    #anhand der Paramlist machen soll !?!
  }

  #update missing rates
  missing_rates <- apply(my_data, 2, function(x) sum(is.na(x)))/nrow(my_data)

  only_one_variable_to_impute <- sum(missing_rates != 0) == 1
  # MS: we have to consider the missing rates in Z too,
  # because we cannot use a random effects variable with NAs.
  missing_rates_Z <- apply(Z, 2, function(x) sum(is.na(x)))/nrow(Z)

  # MS: we have to consider the missing rates in X too,
  # because we cannot use a fixed effects variable with NAs.
  missing_rates_X <- apply(X, 2, function(x) sum(is.na(x)))/nrow(X)

  counter <- 0
  while(any(missing_rates > 0) & counter <= ncol(my_data) + 10){
    #update counter
    counter <- counter + 1

    #find the smallest positive missing rate.
    #If there is not *the one* smallest positive missing rate, sample from those who have the smallest pos. MR.
    tmp_variable <- sample(size = 1,
               names(my_data)[which.min(missing_rates[missing_rates > 0])])

    #get type
    tmp_type <- get_type(my_data[, tmp_variable])

    #generate temporary data set, that contains only the variable to impute and
    #the covariates where no missing values occurs
    tmp_data <- subset(my_data, select = missing_rates == 0 |
                         names(my_data) == tmp_variable)

    tmp_X <- X[, missing_rates_X == 0 & names(X) != tmp_variable, drop = FALSE]
    tmp_Z <- Z[, missing_rates_Z == 0 & names(Z) != tmp_variable, drop = FALSE]


    if(tmp_type == "binary"){

      if(is.null(clID_varname) || !(clID_varname %in% names(tmp_data)) || tmp_variable != multilevel_target_varname){

        imp <- imp_binary(y_variable_name = tmp_variable, data_org_bin_imp = tmp_data,
                        M = M,
                        mn = parameter_list[[tmp_variable]]$mn,
                        max_se = parameter_list[[tmp_variable]]$max_se)
      }else{

        imp <- imp_binary_multi(y_imp_multi = tmp_data[, tmp_variable],
                                X_imp_multi = tmp_X,
                                Z_imp_multi = tmp_Z,
                                clID = tmp_data[, clID_varname],
                                M = M,
                                nitt = 3000,
                                thin = 10,
                                burnin = 1000,
                                allowed_max_value = parameter_list[[tmp_variable]]$allowed_max_value,
                                allowed_max_variable = parameter_list[[tmp_variable]]$allowed_max_variable,
                                allowed_min_value = parameter_list[[tmp_variable]]$allowed_min_value,
                                allowed_min_variable = parameter_list[[tmp_variable]]$allowed_min_variable)



      }
    }

    if(tmp_type == "cont"){

      #MS: take the first "loop" if the data aren't suitable to be imputed by imp.rand()
      #MS: !?!could be improved e.g. by seeting up a boolean.!?! but consider that tmp_data changes
      # and that it could happen,
      #MS: that in the first run a variable is not in tmp_data becaus it still has some NAs
      # but in the second run it is part of tmp_data
      if(is.null(clID_varname) || !(clID_varname %in% names(tmp_data)) || tmp_variable != multilevel_target_varname){

        #now only use the variable with the missing values
        #and the fully observed variables
        #!!!anzahlkinder war character!!!
        #HIER!!!
        imp <- impsyn.cont(y_variable_name = tmp_variable,
                           data.org.cont.imp = tmp_data,
                          cond = parameter_list[[tmp_variable]]$cond,
                          mn = parameter_list[[tmp_variable]]$mn,
                          max_se = parameter_list[[tmp_variable]]$max_se,
                          allowed_max_value = parameter_list[[tmp_variable]]$allowed_max_value,
                          allowed_max_variable = parameter_list[[tmp_variable]]$allowed_max_variable,
                          allowed_min_value = parameter_list[[tmp_variable]]$allowed_min_value,
                          allowed_min_variable = parameter_list[[tmp_variable]]$allowed_min_variable,
                          impsyn = impsyn)
        #!!!brutto einkommen war character!!!

        #If we have hierarchical imputation, use Joerg's Gibbs-sampler
      }else{
        #!!!Die Verwendung des Gibb samplers geht natuerlich nur,
        #wenn mindestens eine random effects variablen samt cluster.id in diesem Schritt verfuegbar sind
        # (d.h. sie tauchen in model_formula auf und sind aktuell vollständig beobachtet und
        # sind deshalb auch in Z_tmp enthalten).

        #Daten mit Gibbs Sampler imputieren
#Gibbs Sampler fuer random coefficitens (also random intercepts und random slopes)
        #?!? und wenn ich eine random effects variable imputieren will,
        #muss ich diese aus Z ausschliessen oder?!?

        imp <- imp_cont_multi(y_imp_multi = tmp_data[, tmp_variable],
                         X_imp_multi = tmp_X,
                         Z_imp_multi = tmp_Z,
                         clID = tmp_data[, clID_varname],
                         M = M,
                         nitt = 3000,
                         thin = 10,
                         burnin = 1000,
                         allowed_max_value = parameter_list[[tmp_variable]]$allowed_max_value,
                         allowed_max_variable = parameter_list[[tmp_variable]]$allowed_max_variable,
                         allowed_min_value = parameter_list[[tmp_variable]]$allowed_min_value,
                         allowed_min_variable = parameter_list[[tmp_variable]]$allowed_min_variable)

       # MS: Wenn ich es sequentiell mache, brauche ich eigentlich nur M = 1???
       #zur sicherheit kann ich auch M = 1 laufen lassen und einfach nur das erste nehmen.

      }
    }

    if(tmp_type == "semicont"){



      if(is.null(clID_varname) || !(clID_varname %in% names(tmp_data)) || tmp_variable != multilevel_target_varname){
        imp <- imp.semi(y_variable_name = tmp_variable, data.org = tmp_data,
                        cond = parameter_list[[tmp_variable]]$cond,
                        mn = parameter_list[[tmp_variable]]$mn,
                        max_se = parameter_list[[tmp_variable]]$max_se,
                        allowed_max_value = parameter_list[[tmp_variable]]$allowed_max_value,
                        allowed_max_variable = parameter_list[[tmp_variable]]$allowed_max_variable,
                        allowed_min_value = parameter_list[[tmp_variable]]$allowed_min_value,
                        allowed_min_variable = parameter_list[[tmp_variable]]$allowed_min_variable,
                        heap = parameter_list[[tmp_variable]]$heap,
                        impsyn = impsyn)
      }else{
        imp <- imp_semicont_multi(y_imp_multi = tmp_data[, tmp_variable],
                              X_imp_multi = tmp_X,
                              Z_imp_multi = tmp_Z,
                              clID = tmp_data[, clID_varname],
                              model_formula = model_formula,
                              heap = 0,
                              M = M,
                              nitt = 3000,
                              thin = 10,
                              burnin = 1000,
                              allowed_max_value = parameter_list[[tmp_variable]]$allowed_max_value,
                              allowed_max_variable = parameter_list[[tmp_variable]]$allowed_max_variable,
                              allowed_min_value = parameter_list[[tmp_variable]]$allowed_min_value,
                              allowed_min_variable = parameter_list[[tmp_variable]]$allowed_min_variable)
      }
    }

    if(tmp_type == "categorical"){
      if(is.null(clID_varname) || !(clID_varname %in% names(tmp_data)) || tmp_variable != multilevel_target_varname){
        imp <- categorical_imp(data.org.categorical.imp = tmp_data, y_variable_name = tmp_variable)

      }else{

        imp <- imp_cat_multi(y_imp_multi = as.factor(tmp_data[, tmp_variable]),
                             X_imp_multi = tmp_X,
                             Z_imp_multi = tmp_Z,
                             clID = tmp_data[, clID_varname],
                             model_formula = model_formula,
                             M = M,
                             nitt = 3000,
                             thin = 10,
                             burnin = 1000,
                             allowed_max_value = parameter_list[[tmp_variable]]$allowed_max_value,
                             allowed_max_variable = parameter_list[[tmp_variable]]$allowed_max_variable,
                             allowed_min_value = parameter_list[[tmp_variable]]$allowed_min_value,
                             allowed_min_variable = parameter_list[[tmp_variable]]$allowed_min_variable)
      }
    }

    if(tmp_type == "ordered_categorical"){
      if(is.null(clID_varname) || !(clID_varname %in% names(tmp_data)) || tmp_variable != multilevel_target_varname){
        imp <- ordered_categorical_imp(data.imp.ordered.categorical = tmp_data, y_variable_name = tmp_variable)
      }else{
        imp <- imp_orderedcat_multi(y_imp_multi = as.factor(tmp_data[, tmp_variable]),
                             X_imp_multi = tmp_X,
                             Z_imp_multi = tmp_Z,
                             clID = tmp_data[, clID_varname],
                             model_formula = model_formula,
                             M = M,
                             nitt = 3000,
                             thin = 10,
                             burnin = 1000,
                             allowed_max_value = parameter_list[[tmp_variable]]$allowed_max_value,
                             allowed_max_variable = parameter_list[[tmp_variable]]$allowed_max_variable,
                             allowed_min_value = parameter_list[[tmp_variable]]$allowed_min_value,
                             allowed_min_variable = parameter_list[[tmp_variable]]$allowed_min_variable)
      }
    }


    if(only_one_variable_to_impute){
      # imputed objects in mice are lists. Each list element is named after a variable in the data set.
      # And each list element contains a nmis times M matrix where the jth column are all nmis imputed data
      # for this variable in the current imputation step.

      imp_hmi <- setNames(list(imp[is.na(my_data[, tmp_variable]), , drop = FALSE]), tmp_variable)

      predictorMatrix <- matrix(0, nrow = ncol(my_data), ncol = ncol(my_data))
      varindex <- colnames(my_data) == tmp_variable
      predictorMatrix[varindex, !varindex]  <- 1

      call <- match.call()
      nmis <- apply(is.na(data), 2, sum)
      state <- list(it = 0, im = 0, co = 0, dep = "", meth = "",
                    log = FALSE)
      loggedEvents <- data.frame(it = 0, im = 0, co = 0, dep = "",
                 meth = "", out = "")

      if (!state$log){ # can state$log ever be TRUE?
        loggedEvents <- NULL
      }
      if (state$log){
        row.names(loggedEvents) <- 1:nrow(loggedEvents)
      }

      if(tmp_type == "categorical" | tmp_type == "ordered_categorical"){

        my_chainMean <- array(dim = M)
        my_chanVar <- array(dim = M)
        for(m in 1:M){
          variable <- (1:nlevels(as.factor(tmp_data[, tmp_variable])))[imp_hmi[[tmp_variable]][,m]]
          my_chainMean[m] <- mean(variable)
          my_chanVar[m] <- var(variable)
        }
      }else{
        my_chainMean <- colMeans(imp_hmi[[tmp_variable]])
        my_chanVar <- apply(imp_hmi[[tmp_variable]], 2, var)
      }

      midsobj <- list(call = call, data = my_data,
                      m = M, nmis = nmis, imp = imp_hmi, method = "multi",
                      predictorMatrix = predictorMatrix, visitSequence = varindex,
                      form = form, post = post, seed = NA, iteration = NA,
                      lastSeedValue = .Random.seed,
                      chainMean = my_chainMean,
                      chainVar = my_chanVar,
                      loggedEvents = loggedEvents)
#!!!NEED TO ADD vcov() method to midsobj!!!
      oldClass(midsobj) <- "mids"
      return(midsobj)
    }


    #MS: update my_data with imputed variable
    my_data[, tmp_variable] <- imp

    #MS: update Z, because it could be the case that a variable contained in Z
    #MS: was udated by imputation
    if(tmp_variable %in% names(Z)) Z[, tmp_variable] <- imp
    #MS: update X, because it could be the case that a variable contained in Z
    #MS: was udated by imputation
    if(tmp_variable %in% names(X)) X[, tmp_variable] <- imp

    #update missing rates
    missing_rates <- apply(my_data, 2, function(x) sum(is.na(x)))/nrow(my_data)
    missing_rates_Z <- apply(Z, 2, function(x) sum(is.na(x)))/nrow(Z)
    missing_rates_X <- apply(X, 2, function(x) sum(is.na(x)))/nrow(X)

  }

  ntries <- 1000#Number of tries to reach convergence

  convergence <- FALSE
  number.produced.datasets <- 0
  counter <- 0
  #MS: laufe solange durch, bis M datens?tze produziert wurden, oder die maximale Anzahl an Iterationen erreicht wurde.
  #MS: um M Datens?tze zu produzieren, m?ssen wir warten, bis bei allen Variablen der burn in erreicht ist.
  #MS: Danach m?ssen wir gem?? der ermittelten notwendigen Schrittweite s jeden s-te Datensatz ausgegeben,
  #MS: bis insgesamt M Datens?tze vorliegen.
  data_of_all_iterations <- list()
  while(number.produced.datasets < M & counter < ntries){

    counter <- counter + 1

    for(l2 in incomplete.variables){


      #get type
      #Hinweis: get.type wird eigentlich nicht mehr gebraucht, es sei denn
      #eine semicontinuierliche variable wird nur noch als cont eingestuft.
      tmp_type <- get.type(my_data[, l2])#!!!anzahlkinder ist aus irgendeinem Grund
      #ein character, weshalb get.type() nicht mehr einen typ ermitteln kann!!!

      #default value for imp
      imp <- my_data[, l2]#?!? Sollte aber eigentlich nicht gebraucht werden
      #sofern eine der if-Schleifen eintritt (was schon sehr wichtig w?re)
      #diese Default-Wert-Zuweisung ist also als eine Art Sicherheitszeile.

      #make currrent variable which is about to be imputed again containing
      #the original missing values
      my_data[, l2] <- data[, l2]
      #Hier w?re es also kein Problem, wenn man tmp_data so stricken w?rde, dass
      #die gefilterten Fragen faktorisiert werden.
      #Die R-Variablen am besten in data, data.work, data.ret etc. umbenennen,
      #so dass klar ist: was ist orginal, welche daten werden nur f?r die Bearbeitung
      #manipuliert und welche werden am Ende ausgegeben?



      if(tmp_type == "binary"){

        imp <- imp_binary(y_variable_name = l2, data_org_bin_imp = my_data,
                          cond = parameter_list[[l2]]$cond,
                          mn = parameter_list[[l2]]$mn,
                          max_se = parameter_list[[l2]]$max_se,
                          impsyn = impsyn)
      }

      if(tmp_type == "cont"){


        #If we have hierarchical imputation, use J?rg's Gibbs-sampler


        if(is.null(clID_varname) || !(clID_varname %in% names(my_data)) || l2 != target_varname){


          #now only use the variable with the missing values
          #and the fully observed variables

          imp <- impsyn.cont(y_variable_name = l2, data.org.cont.imp = my_data,
                            cond = parameter_list[[l2]]$cond,
                            mn = parameter_list[[l2]]$mn,
                            max_se = parameter_list[[l2]]$max_se,
                            allowed_max_value = parameter_list[[l2]]$allowed_max_value,
                            allowed_max_variable = parameter_list[[l2]]$allowed_max_variable,
                            allowed_min_value = parameter_list[[l2]]$allowed_min_value,
                            allowed_min_variable = parameter_list[[l2]]$allowed_min_variable,
                            impsyn = impsyn)

        }else{

  #TODO: MAKE SHURE, THAT Z ONLY CONTAINS VARIABLES THAT ARE ASSUMED TO HAVE A RANDOM EFFECT!!!
		      #Gibbs Sampler f?r random coefficitens (also random intercepts und random slopes)
	        imp <- impsyn.cont(y_variable_name = l2, data.org.cont.imp = my_data,
	                           cond = parameter_list[[l2]]$cond,
	                           mn = parameter_list[[l2]]$mn,
	                           max_se = parameter_list[[l2]]$max_se,
	                           allowed_max_value = parameter_list[[l2]]$allowed_max_value,
	                           allowed_max_variable = parameter_list[[l2]]$allowed_max_variable,
	                           allowed_min_value = parameter_list[[l2]]$allowed_min_value,
	                           allowed_min_variable = parameter_list[[l2]]$allowed_min_variable,
	                           impsyn = impsyn)

      #MS: !!!FUERS SCHNELLERE TESTEN, WURDE HIER DER GIBBS SAMPLER DURCH DIE NORMALE CONT-IMPUTATION ERSETZT!!!

            #imp.rand(y_variable_name = l2, data.imp.rand = my_data,
	         #       Z.imp.rand = tmp.Z,
	#                cl.id =  subset(my_data, select = clID_varname),
	#                n.iter = 100, m = 10, n.chains = 3, burn.in = 1/3)
        }
      }

      if(tmp_type == "semicont"){

        imp <- imp.semi(y_variable_name = l2, data.org = my_data,
                        cond = parameter_list[[l2]]$cond,
                        mn = parameter_list[[l2]]$mn,
                        max_se = parameter_list[[l2]]$max_se,
                        allowed_max_value = parameter_list[[l2]]$allowed_max_value,
                        allowed_max_variable = parameter_list[[l2]]$allowed_max_variable,
                        allowed_min_value = parameter_list[[l2]]$allowed_min_value,
                        allowed_min_variable = parameter_list[[l2]]$allowed_min_variable,
                        heap = parameter_list[[l2]]$heap,
                        impsyn = impsyn)
      }

      if(tmp_type == "categorical"){
        imp <- categorical_imp(data.org.categorical.imp = my_data, y_variable_name = l2)
      }

      if(tmp_type == "ordered_categorical"){
        imp <- ordered_categorical_imp(data.imp.ordered.categorical = my_data, y_variable_name = l2)
      }


      # MS: it could be the case that the user wants to impute only one variable.
      # In this case I just can return the M data sets

      my_data[, l2] <- imp


      #get theta (up to now: the mean) of the impuation sequence
      #?what to use for the categorical values?
      #I excluded the non imputed values to compare only values that *could* have
      #changed after the last imputation

      #since I have no idea what to measure for a categorical variable
      #"whats the mean of c("hans", "berta", "berta", "hans")?
      #use the mode instead? what would be the mode here?
      #what would be the mode of a continous variable?
      #I decided to estimate the influence of the current variable on a continous
      #one in a simple linear model. and i only look onto the intercept because
      #this exist for every variable. It doesnt have to be meaningfull.
      #It only has to remain equal for equal data.
      #and small changes in the data should lead to only small in the estimate
      #(the same for "big").

      #MS: update Z, because it could be the case that a variable contained in Z
      #MS: was udated by imputation
      if(l2 %in% names(Z)) Z[, l2] <- imp


      #MS: TODO: AUF BURN-IN TESTEN (NACH GEWCKE oder WELCH) HIER WEITER!!!
      #MS: DATENSATZ SPEICHERN UND AUF KONVERGENZ TESTEN

      #MS: TODO: USE MORE ELABORATED COEFFICIENTS (SEE HEAPING CODE)!!!
      # mean(imp[is.na(data[, tmp_variable]), ])
    }#MS: end of l2 in incomplete.variables loop

    #MS: aktuellen datensatz abspeichern

    data_of_all_iterations[[counter]] <- my_data

    #MS: check for convergence, but only if convergence never has been TRUE
    #MS: Der Konvergenz-Check und die Burnin-Ermittlung benoetigen die ValuesToMonitor.
    #MS: Die werden auch nach Konvergenz-Check und Burnin-Ermmitlungen gebraucht um den lag zu ermitteln,
    #MS: weshalb sie immer wieder neu ermittelt werden - aber der Burnin muss weg geschmissen werden.

    #check convergence
    convergence.tmp <- convergence_diagnostics(daten.convergence_diagnostics = data_of_all_iterations[-(1:burnin)],
                                               interestingvariables = incomplete.variables,
                                               nonresponse_matrix = nonresponse_matrix,
                                               parameter_list = parameter_list,
                                               criterion = "HeidelbergerWelch",
                                               heap = 0,
                                               plotit = FALSE)

    #MS: Bei der Heidelberger and Welch Diagnostic, wird eine Liste ausgegeben
    if(is.list(convergence.tmp))  ValuesToMonitor <- convergence.tmp$ValuesToMonitor

    if(!convergence){
      #MS: Bei der Heidelberger and Welch Diagnostic, wird eine Liste ausgegeben
      if(is.list(convergence.tmp)){
        burnin <- convergence.tmp$burnin
        convergence <- all(convergence.tmp$converged_values)

      }else{
        convergence <- all(convergence.tmp)
      }

    }else{ #MS: Fall, dass Konvergenz eingetreten ist.


      #MS: check lag after discarding the burn in.
      lagg <- check.lag(estimates = (ValuesToMonitor[, 1]),
                        #MS: !!!jitter() entfernen und auf jede Spalte von valuesToMonitor anwenden.!!!
                        burnin = burnin,
                        plotit = FALSE)
      #MS: check.lag() braucht ValuesToMonitor, und das wird momentan nur ausgegeben,
      #MS: wenn bei convergence_diagnostics() der Parameter criterion = "HeidelbergerWelch" uebergeben wird.

      #MS: ermitteln, welche Datensaetze verwendet werden koennen
      verwertbar <- seq(from = burnin, to = length(data_of_all_iterations), by = lagg)
      #?!? lagg updaten???
      #MS: lagg is intentionally mispelled to avoid confusion with the function lag()

      #MS: die Anzahl der verwertbaren Datensaetze ergibt sich dann aus der Laenge der Indizes.
      number.produced.datasets <- length(verwertbar)
    }
    #MS: !!! TODO: NACH KONVERGENZ BURN IN ERMITTELN, SCHRITTWEITE ERMITETELN,
    #MS: UND DANN NOCH ENTSPRECHEND OFT DIE ITERATIONEN LAUFEN LASSEN!!!

  }#MS: Ende number.produced.datasets < M & counter < ntries while-Schleife
#http://en.wikipedia.org/wiki/Trace_(linear_algebra)
  #MS: jetzt die ausgewaehlten Datensaetze ausgeben.
  return(data_of_all_iterations[verwertbar])

}



