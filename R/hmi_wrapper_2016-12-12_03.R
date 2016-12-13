
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
#' @param maxit An integer defining the number of times the imputation cycle
#' (imputing x_1|x_{-1} than x_2|x_{-2}, ... x_p|x_{-p}) shall be repeated.
#' The task of cheching convergence is left to the user, by evaluating the chainMean and chainVar!
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
                    maxit = 5){

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

  #get missing rates
  missing_rates <- apply(my_data, 2, function(x) sum(is.na(x)))/nrow(my_data)

  only_one_variable_to_impute <- sum(missing_rates != 0) == 1

  #MS: it's not necessary to divide by nrow(my_data)

  #get variables with missings
  incomplete_variables <- names(my_data)[missing_rates > 0]

  #MS: Nonresponse-Matrix aufbauen: (wird später für die Convergenz Diagnostic gebraucht)
  nonresponse_matrix <- is.na(my_data)


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
    multilevel_target_varname_full <- as.character(model_formula)[2]

    #check if the formula contains any functions like "log(y)" in log(y) ~ x.
    multilevel_target_varname <- gsub(".*\\((.*)\\).*", "\\1", multilevel_target_varname_full)

    if(multilevel_target_varname != multilevel_target_varname_full){
      warning(paste("For the tool to work saver we suggest to do the transformation -->", multilevel_target_varname_full,
                    "<-- before running the imputation."))
    }

    # check if this multilevel_target_varname realy occurs in the dataest
    if(!multilevel_target_varname %in% names(my_data)){
      stop(paste(multilevel_target_varname, "not found in dataset. Please check your model_formula and data."))
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



  }else{ # if model_formula is NULL
    # in this case impute every variable based on the others in a single level framework
	  clID_varname <- NULL
	  fixedeffects_varname <- colnames(my_data)
	  randomeffects_varname <- character()
  }


  #In order to have a rectangular data.set (without NAs), do the sample_imp:
  for(j in 1:ncol(my_data)){
    my_data[, j] <- sample_imp(my_data[, j])#apply would be smarter, but it returns a matrix.
  }


  #MS: now safe the variables used for modeling fixed effects
  X <- my_data[, fixedeffects_varname, drop = FALSE]

  #MS: now safe the variables used for modeling random effects
  Z <- my_data[, randomeffects_varname, drop = FALSE]

  #Preparation for setting up a mids-object
  alldata <- list()#this object saves the M datasets.
  my_chainMean <- array(dim = c(length(incomplete_variables), maxit, M))
  dimnames(my_chainMean) <-  list(incomplete_variables,
                                  as.character(1:maxit),
                                  paste("Chain", 1:M))
  my_chainVar <- my_chainMean
  for(i in 1:M){
    for(l1 in 1:maxit){

      for(l2 in incomplete_variables){#do the imputation cycle
        #find the smallest positive missing rate.
        #If there is not *the one* smallest positive missing rate, sample from those who have the smallest pos. MR.
        tmp_variable <- l2

        #get type
        tmp_type <- get_type(my_data[, tmp_variable])

        #generate temporary data set, that contains only the variable to impute and
        #the covariates where no missing values occurs
        #tmp_data <- my_data

        #make currrent variable which is about to be imputed again containing
        #the original missing values
        my_data[, l2] <- data[, l2]

        tmp_X <- X[, names(X) != tmp_variable, drop = FALSE]
        tmp_Z <- Z[, names(Z) != tmp_variable, drop = FALSE]


        if(tmp_type == "binary"){

          if(is.null(clID_varname) || !(clID_varname %in% names(my_data)) || tmp_variable != multilevel_target_varname){

            imp <- imp_binary(y_imp_multi = my_data[, tmp_variable],
                             X_imp_multi = tmp_X,
                             M = 1,
                             allowed_max_value = parameter_list[[tmp_variable]]$allowed_max_value,
                             allowed_max_variable = parameter_list[[tmp_variable]]$allowed_max_variable,
                             allowed_min_value = parameter_list[[tmp_variable]]$allowed_min_value,
                             allowed_min_variable = parameter_list[[tmp_variable]]$allowed_min_variable)

          }else{

            imp <- imp_binary_multi(y_imp_multi = my_data[, tmp_variable],
                                    X_imp_multi = tmp_X,
                                    Z_imp_multi = tmp_Z,
                                    clID = my_data[, clID_varname],
                                    M = 1,
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

          if(is.null(clID_varname) || !(clID_varname %in% names(my_data)) || tmp_variable != multilevel_target_varname){


            imp <- imp_cont(y_imp_multi = my_data[, tmp_variable],
                                  X_imp_multi = tmp_X,
                                  M = 1,
                                  allowed_max_value = parameter_list[[tmp_variable]]$allowed_max_value,
                                  allowed_max_variable = parameter_list[[tmp_variable]]$allowed_max_variable,
                                  allowed_min_value = parameter_list[[tmp_variable]]$allowed_min_value,
                                  allowed_min_variable = parameter_list[[tmp_variable]]$allowed_min_variable)

          }else{


            imp <- imp_cont_multi(y_imp_multi = my_data[, tmp_variable],
                             X_imp_multi = tmp_X,
                             Z_imp_multi = tmp_Z,
                             clID = my_data[, clID_varname],
                             M = 1,
                             nitt = 3000,
                             thin = 10,
                             burnin = 1000,
                             allowed_max_value = parameter_list[[tmp_variable]]$allowed_max_value,
                             allowed_max_variable = parameter_list[[tmp_variable]]$allowed_max_variable,
                             allowed_min_value = parameter_list[[tmp_variable]]$allowed_min_value,
                             allowed_min_variable = parameter_list[[tmp_variable]]$allowed_min_variable)

          }
        }

        if(tmp_type == "semicont"){



          if(is.null(clID_varname) || !(clID_varname %in% names(my_data)) || tmp_variable != multilevel_target_varname){
            imp <- imp_semicont(y_imp_multi = my_data[, tmp_variable],
                                      X_imp_multi = tmp_X,
                                      heap = 0,
                                      M = 1,
                                      allowed_max_value = parameter_list[[tmp_variable]]$allowed_max_value,
                                      allowed_max_variable = parameter_list[[tmp_variable]]$allowed_max_variable,
                                      allowed_min_value = parameter_list[[tmp_variable]]$allowed_min_value,
                                      allowed_min_variable = parameter_list[[tmp_variable]]$allowed_min_variable)
          }else{
            imp <- imp_semicont_multi(y_imp_multi = my_data[, tmp_variable],
                                  X_imp_multi = tmp_X,
                                  Z_imp_multi = tmp_Z,
                                  clID = my_data[, clID_varname],
                                  model_formula = model_formula,
                                  heap = 0,
                                  M = 1,
                                  nitt = 3000,
                                  thin = 10,
                                  burnin = 1000,
                                  allowed_max_value = parameter_list[[tmp_variable]]$allowed_max_value,
                                  allowed_max_variable = parameter_list[[tmp_variable]]$allowed_max_variable,
                                  allowed_min_value = parameter_list[[tmp_variable]]$allowed_min_value,
                                  allowed_min_variable = parameter_list[[tmp_variable]]$allowed_min_variable)
          }
        }


        if(tmp_type == "roundedcont"){
          if(is.null(clID_varname) || !(clID_varname %in% names(my_data)) || tmp_variable != multilevel_target_varname){

            imp <- imp_roundedcont(y_imp_multi = my_data[, tmp_variable],
                                         X_imp_multi = tmp_X,
                                         intercept_varname = NULL,
                                         M = 1,
                                         allowed_max_value = parameter_list[[tmp_variable]]$allowed_max_value,
                                         allowed_max_variable = parameter_list[[tmp_variable]]$allowed_max_variable,
                                         allowed_min_value = parameter_list[[tmp_variable]]$allowed_min_value,
                                         allowed_min_variable = parameter_list[[tmp_variable]]$allowed_min_variable)



          }else{
            #print("If you want to impute a (rounded) income variable, consider to take the logarithm of it and to standardize your covariates.")
            #???CAN WE DO IT IN A MULTILEVEL STYLE???
            imp <- imp_roundedcont_multi(y_imp_multi = my_data[, tmp_variable],
                                         X_imp_multi = tmp_X,
                                         Z_imp_multi = tmp_Z,
                                         clID = my_data[, clID_varname],
                                         intercept_varname = intercept_varname,
                                         M = 1,
                                         allowed_max_value = parameter_list[[tmp_variable]]$allowed_max_value,
                                         allowed_max_variable = parameter_list[[tmp_variable]]$allowed_max_variable,
                                         allowed_min_value = parameter_list[[tmp_variable]]$allowed_min_value,
                                         allowed_min_variable = parameter_list[[tmp_variable]]$allowed_min_variable)


          }
        }


        if(tmp_type == "categorical"){
          if(is.null(clID_varname) || !(clID_varname %in% names(my_data)) || tmp_variable != multilevel_target_varname){

            imp <- imp_cat(data.org.categorical.imp = my_data, y.variable.name = tmp_variable)

          }else{

            imp <- imp_cat_multi(y_imp_multi = as.factor(my_data[, tmp_variable]),
                                 X_imp_multi = tmp_X,
                                 Z_imp_multi = tmp_Z,
                                 clID = my_data[, clID_varname],
                                 model_formula = model_formula,
                                 M = 1,
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
          if(is.null(clID_varname) || !(clID_varname %in% names(my_data)) || tmp_variable != multilevel_target_varname){
            imp <- imp_orderedcat(data.imp.ordered.categorical = my_data, y.variable.name = tmp_variable,
                                  M = 1)
          }else{
            imp <- imp_orderedcat_multi(y_imp_multi = as.factor(my_data[, tmp_variable]),
                                 X_imp_multi = tmp_X,
                                 Z_imp_multi = tmp_Z,
                                 clID = my_data[, clID_varname],
                                 model_formula = model_formula,
                                 M = 1,
                                 nitt = 3000,
                                 thin = 10,
                                 burnin = 1000,
                                 allowed_max_value = parameter_list[[tmp_variable]]$allowed_max_value,
                                 allowed_max_variable = parameter_list[[tmp_variable]]$allowed_max_variable,
                                 allowed_min_value = parameter_list[[tmp_variable]]$allowed_min_value,
                                 allowed_min_variable = parameter_list[[tmp_variable]]$allowed_min_variable)
          }
        }



        #MS: update my_data with imputed variable
        my_data[, tmp_variable] <- imp

        #MS: update Z, because it could be the case that a variable contained in Z
        #MS: was udated by imputation
        if(tmp_variable %in% names(Z)) Z[, tmp_variable] <- imp
        #MS: update X, because it could be the case that a variable contained in Z
        #MS: was udated by imputation
        if(tmp_variable %in% names(X)) X[, tmp_variable] <- imp

        to_evaluate <- imp[is.na(data[, tmp_variable]), ]

        if(tmp_type == "categorical" | tmp_type == "ordered_categorical"){

          to_evaluate_numerical <- (1:nlevels(as.factor(my_data[, tmp_variable])))[to_evaluate]

        }else{
          to_evaluate_numerical <- to_evaluate
        }

        my_chainMean[l2, l1, i]  <- mean(to_evaluate_numerical)
        my_chainVar[l2, l1, i]  <- var(to_evaluate_numerical)
      }

    }
    alldata[[i]] <- my_data
  }

  ############################        SET UP MIDS OBJECT            ####################

#???HOW DOES THE MIDS OBJECT HAS TO LOOK IF THERE ARE MORE THAN ONE VARIABLE TO IMPUTE???
  if(only_one_variable_to_impute){
    # imputed objects in mice are lists. Each list element is named after a variable in the data set.
    # And each list element contains a nmis times M matrix where the jth column are all nmis imputed data
    # for this variable in the current imputation step.

    imp_hmi <- setNames(list(imp[is.na(my_data[, tmp_variable]), , drop = FALSE]), tmp_variable)

    #predictorMatrix <- matrix(0, nrow = ncol(my_data), ncol = ncol(my_data))
    #varindex <- colnames(my_data) == tmp_variable
    #predictorMatrix[varindex, !varindex]  <- 1

    call <- match.call()

    state <- list(it = 0, im = 0, co = 0, dep = "", meth = "",
                  log = FALSE)


    if (!state$log){ # can state$log ever be TRUE?
      loggedEvents <- NULL
    }
    if (state$log){
      row.names(loggedEvents) <- 1:nrow(loggedEvents)
    }



  }else{ # if there is more than one variable to return



    imp_hmi <- setNames(vector("list", ncol(data)), colnames(data))
    for(l2 in incomplete_variables){
      n_mis_j <- sum(is.na(data[, l2]))
      imp_hmi[[l2]] <- as.data.frame(matrix(NA, nrow = n_mis_j, ncol = M, dimnames = list(NULL, 1:M)))
      for(i in 1:M){
        imputed_values <- alldata[[i]][l2][is.na(data[, l2]), , drop = FALSE]
        imp_hmi[[l2]][, i] <- imputed_values
        rownames(imp_hmi[[l2]]) <- rownames(imputed_values)
      }
    }

  }

  nmis <- apply(is.na(data), 2, sum)
  method <- rep("multi", times = ncol(data)) #!!! DAS  MUSS ICH NOCH AUSBESSERN JE NACH DEM OB DIE METHODE WIRKLICH MULTI WAR ODER NICHT
  #UND AUCH DANACH OB ES CONT; SEMICONT; CAT ETC WAR!!!
  names(method) <- colnames(data)

  predictorMatrix <- matrix(0, nrow = ncol(my_data), ncol = ncol(my_data))
  rownames(predictorMatrix) <- colnames(data)
  colnames(predictorMatrix) <- colnames(data)

  predictorMatrix[incomplete_variables, ] <- 1
  diag(predictorMatrix) <- 0

  visitSequence <- (1:ncol(data))[apply(is.na(data), 2, any)]
  names(visitSequence) <- colnames(data)[visitSequence]

  loggedEvents <- data.frame(it = 0, im = 0, co = 0, dep = "",
                             meth = "", out = "")

  #Not supported yet
  form <- vector("character", length = ncol(data))
  post <- vector("character", length = ncol(data))

  midsobj <- list(call = call, data = data,
                  m = M, nmis = nmis, imp = imp_hmi, method = method,
                  predictorMatrix = predictorMatrix, visitSequence = visitSequence,
                  form = form, post = post, seed = NA, iteration = maxit,
                  lastSeedValue = .Random.seed,
                  chainMean = my_chainMean,
                  chainVar = my_chainVar,
                  loggedEvents = loggedEvents)
  #!!!NEED TO ADD vcov() method to midsobj!!!
  oldClass(midsobj) <- "mids"
  return(midsobj)

}



