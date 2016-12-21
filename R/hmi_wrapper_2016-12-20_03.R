
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
#' hmi_imp <- hmi(data = test, model_formula = my.formula)
#' hmifit <- hmi_pool(data = hmi_imp, analysis = my_analysis)
#' pool(with(data = hmi_imp, expr = lmer(Reaction ~ Days + (1 + Days | Subject))))

hmi <- function(data, model_formula = NULL,
                    M = 10,
                    nitt = 3000,
                    thin = 100,
                    burnin = 1000,
                    maxit = 5){

  my_data <- data #wird spaeter ausgegeben. Entspricht von der Form dem
  #original Datensatz, nur, dass alle fehlenden Werte (mit Ausnahme der design
  #bedingt fehlenden) durch imputierte Werte ersetzt wurden.

  #des weitern gibt es noch einen tmp_data Datensatz

  if(is.matrix(my_data)) stop("We need your data in a data.frame (and not a matrix).")


  if(!is.data.frame(my_data)) stop("We need your data in a data.frame.")
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

  # # # # # # # # # # #  get the formula elements (fe) and do consistency checks  # # # # # # # # # # # # # # # #
  constant_variables <- apply(my_data, 2, function(x) length(unique(x)) == 1)

  if(sum(constant_variables) >= 2){
    stop(paste("Your data has two or more constant variables: ", names(my_data)[constant_variables],
               ". One needs to be dropped to avoid multicolinearity.", sep = ""))
  }

  variable_names_in_data <- colnames(my_data)
  fe <- extract_varnames(model_formula = model_formula, constant_variables = constant_variables,
                   variable_names_in_data = variable_names_in_data)


  #check if this variable name realy occurs in the dataset
  if(all(fe$fixedeffects_varname != "") & !all(fe$fixedeffects_varname %in% names(my_data))){
    writeLines(paste("We didn't find",
                  paste(fe$fixedeffects_varname[!fe$fixedeffects_varname %in% names(my_data)],
                        collapse = " and "),
                  "in your data. How do you want to proceed: \n
[c]ontinue with ignoring the model_formula an running a single level imputation
or [e]xiting the imputation?"))

    proceed <- readline("Type 'c' or 'e' into the console and press [enter]: ")
    if(proceed == "e") return(NULL)
    #if the model_formula is ignored every extracted value of the model_formula is ignored...
    model_formula <- NULL
    fe[1:length(fe)] <- ""
    #... and every variable is used as a fixed effects variable, because in this case a single level imputation
    #on every variable is run
    fe$fixedeffects_varname <- colnames(my_data)
  }

  #check if this variable name realy occurs in the dataset
  if(all(fe$randomeffects_varname != "") & !all(fe$randomeffects_varname %in% names(my_data))){

    writeLines(paste("We didn't find",
                     paste(fe$randomeffects_varname[!fe$randomeffects_varname %in% names(my_data)],
                           collapse = " and "),
                     "in your data. How do you want to proceed: \n
                     [c]ontinue with ignoring the model_formula an running a single level imputation
                     or [e]xiting the imputation?"))

    proceed <- readline("Type 'c' or 'e' into the console and press [enter]: ")


    if(proceed == "e") return(NULL)
    model_formula <- NULL
    fe[1:length(fe)] <- ""
    fe$fixedeffects_varname <- colnames(my_data)
  }



  if(all(c("0", "1") %in% fe$randomeffects_varname)){
    warning("Your model_formula includes 0 for no random intercept and 1 for a random intercept at the same time.
            Please decide whether you want a random intercept or not.")
  }


  #it could be the case that the user specifies several cluster ids.
  #As the tool doesn't support this yet, we stop the imputation routine.
  if(length(fe$clID_varname) >= 2){
    writeLines("The package only supports one cluster ID in the model_formula. \n
How do you want to proceed: \n
                     [c]ontinue with ignoring the model_formula an running a single level imputation
                     or [e]xiting the imputation?")

    proceed <- readline("Type 'c' or 'e' into the console and press [enter]: ")
    if(proceed == "e") return(NULL)
    model_formula <- NULL
    fe[1:length(fe)] <- ""
    fe$fixedeffects_varname <- colnames(my_data)
  }

  #check whether the formula matches the data
  if(fe$clID_varname != "" & !(fe$clID_varname %in% names(my_data))){
    writeLines(paste("We didn't find", fe$clID_varname,
                     "in your data. How do you want to proceed: \n
                     [c]ontinue with ignoring the model_formula an running a single level imputation
                     or [e]xiting the imputation?"))

    proceed <- readline("Type 'c' or 'e' into the console and press [enter]: ")
    if(proceed == "e") return(NULL)
    model_formula <- NULL
    fe[1:length(fe)] <- ""
    fe$fixedeffects_varname <- colnames(my_data)
  }

  #check whether the constant variable in the dataset are labbeld "1" if this is the label of the
  #intercept variable in the model formula.
  if(any(constant_variables)){
    if((any(fe$fixedeffects_varname == "1") &
        variable_names_in_data[constant_variables] != "1")|
       (any(fe$randomeffects_varname == "1") &
        variable_names_in_data[constant_variables] != "1")){
      warning(paste("Your model_formula includes '1' for an intercept.
                          The constant variable in your data is named ", variable_names_in_data[constant_variables], ".
                          Consider naming the variable also '1' or 'Intercept'.", sep = ""))
    }
  }


  #make shure that the random effects variable exist, and has the value 1
  if(fe$intercept_varname != ""){
    print(paste("We interprete", fe$intercept_varname,
                "as the intercept variable and set its value to 1."))
    my_data[, fe$intercept_varname] <- 1
  }

  if(fe$clID_varname != ""){
    print(paste("We interprete", fe$clID_varname,
                "as the cluster indicator and treat it as a factor."))
    my_data[, fe$clID_varname] <- as.factor(my_data[, fe$clID_varname])
  }



  ###########################  MERGE SMALL CLUSTERS ####################
  if(fe$clID_varname != ""){

    tab_1 <- table(my_data[, fe$clID_varname])

    # set a minumim size for the clusters
    mn <- 6

    #merge small clusters to a large big cluster
    if(any(tab_1  < mn)){
      writeLines(paste(sum(tab_1 < mn), "clusters have less than", mn, "observations.",
                       "How do you want to proceed: \n
                       [c]ontinue with merging these clusters with others
                       or [e]xiting the imputation?"))

      proceed <- readline("Type 'c' or 'e' into the console and press [enter]: ")
      if(proceed == "e") return(NULL)

      # The procedure is the following:
      # 1. determine the smallest and the second smallest cluster.
      # 2. take the observations of the smallest cluster and combine it with observations of the second smallest cluster.
      #  this is be done by giving both clusters the name "[name of second smallest cluster] forced merge".
      # 3. repeat with 1 until there is no invalid cluster left
      while(any(tab_1 < mn)){

        #step 1: determine the smallest and the second smallest cluster.
        smallest_cluster <- tab_1[tab_1 == sort(as.numeric(tab_1))[1]][1]
        tab_1_without_smallest_cluster <- tab_1[names(tab_1) != names(smallest_cluster)]
        second_smallest_cluster <- tab_1_without_smallest_cluster[tab_1_without_smallest_cluster == sort(as.numeric(tab_1_without_smallest_cluster))[1]][1]

        #step 2: new cluster name
        new_name <-  paste(names(second_smallest_cluster), "forced_merge", sep = "_")

        levels(my_data[, fe$clID_varname])[levels(my_data[, fe$clID_varname]) == names(smallest_cluster)] <- new_name
        levels(my_data[, fe$clID_varname])[levels(my_data[, fe$clID_varname]) == names(second_smallest_cluster)] <- new_name
        #step 3: repeat.

        tab_1 <- table(my_data[, fe$clID_varname])

      }

    }

  }

  ########################################   START THE IMPUTATION #####################

  #In order to have a rectangular data.set (without NAs), do the sample_imp:
  types <- array(dim = ncol(my_data))#get the variable types

  for(j in 1:ncol(my_data)){
    my_data[, j] <- sample_imp(my_data[, j])#apply would be smarter, but it returns a matrix.
    types[j] <- get_type(my_data[, j])
  }

  categorical <- types == "categorical"
  if(any(categorical)){
    for(catcount in 1:sum(categorical)){
      print(paste("We interprete", names(my_data)[categorical][catcount],
                  "as the categorical variable and force it to be a factor."))
      my_data[, names(my_data)[categorical][catcount]] <-
        as.factor(my_data[, names(my_data)[categorical][catcount]])

    }
  }

  #MS: now safe the variables used for modeling fixed effects

  X <- my_data[, fe$fixedeffects_varname[fe$fixedeffects_varname %in% names(my_data)], drop = FALSE]

  #MS: now safe the variables used for modeling random effects
  Z <- my_data[, fe$randomeffects_varname[fe$randomeffects_varname %in% names(my_data)] , drop = FALSE]





  ####################   Preparation for setting up a mids-object
  alldata <- list()#this object saves the M datasets.
  my_chainMean <- array(dim = c(length(incomplete_variables), maxit, M))
  dimnames(my_chainMean) <-  list(incomplete_variables,
                                  as.character(1:maxit),
                                  paste("Chain", 1:M))

  my_chainVar <- my_chainMean
  for(i in 1:M){
    for(l1 in 1:maxit){

      for(l2 in incomplete_variables){#do the imputation cycle
        print(Sys.time())
        print(paste("Imputing variable ", l2, " at iteration ", l1,
                    " for imputation ", i, " (out of ", M, ")", sep = ""))
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

        #check number of observed values in each cluster.

        if(fe$clID_varname != ""){
          tab_1 <- table(my_data[!is.na(my_data[, tmp_variable]), fe$clID_varname])

          #merge small clusters to a large big cluster
          while(any(tab_1 < mn)){

            #step 1: determine the smallest and the second smallest cluster.
            smallest_cluster <- tab_1[tab_1 == sort(as.numeric(tab_1))[1]][1]
            tab_1_without_smallest_cluster <- tab_1[names(tab_1) != names(smallest_cluster)]
            second_smallest_cluster <- tab_1_without_smallest_cluster[tab_1_without_smallest_cluster == sort(as.numeric(tab_1_without_smallest_cluster))[1]][1]

            #step 2: new cluster name
            new_name <-  paste(names(second_smallest_cluster), "forced_merge", sep = "_")

            levels(my_data[, fe$clID_varname])[levels(my_data[, fe$clID_varname]) == names(smallest_cluster)] <- new_name
            levels(my_data[, fe$clID_varname])[levels(my_data[, fe$clID_varname]) == names(second_smallest_cluster)] <- new_name
            #step 3: repeat.

            tab_1 <- table(my_data[!is.na(my_data[, tmp_variable]), fe$clID_varname])

          }

        }



        if(tmp_type == "binary"){

          if(is.null(fe$clID_varname) || !(fe$clID_varname %in% names(my_data)) || tmp_variable != fe$target_varname){

            imp <- imp_binary(y_imp_multi = my_data[, tmp_variable],
                             X_imp_multi = tmp_X,
                             M = 1)

          }else{

            imp <- imp_binary_multi(y_imp_multi = my_data[, tmp_variable],
                                    X_imp_multi = tmp_X,
                                    Z_imp_multi = tmp_Z,
                                    clID = my_data[, fe$clID_varname],
                                    M = 1,
                                    nitt = nitt,
                                    thin = thin,
                                    burnin =  burnin)



          }
        }

        if(tmp_type == "cont"){

          if(is.null(fe$clID_varname) || !(fe$clID_varname %in% names(my_data)) || tmp_variable != fe$target_varname){


            imp <- imp_cont(y_imp_multi = my_data[, tmp_variable],
                                  X_imp_multi = tmp_X,
                                  M = 1)

          }else{


            imp <- imp_cont_multi(y_imp_multi = my_data[, tmp_variable],
                             X_imp_multi = tmp_X,
                             Z_imp_multi = tmp_Z,
                             clID = my_data[, fe$clID_varname],
                             M = 1,
                             nitt = 3000,
                             thin = 10,
                             burnin = 1000)

          }
        }

        if(tmp_type == "semicont"){



          if(is.null(fe$clID_varname) || !(fe$clID_varname %in% names(my_data)) || tmp_variable != fe$target_varname){
            imp <- imp_semicont(y_imp_multi = my_data[, tmp_variable],
                                      X_imp_multi = tmp_X,
                                      heap = 0,
                                      M = 1)
          }else{
            imp <- imp_semicont_multi(y_imp_multi = my_data[, tmp_variable],
                                  X_imp_multi = tmp_X,
                                  Z_imp_multi = tmp_Z,
                                  clID = my_data[, fe$clID_varname],
                                  model_formula = model_formula,
                                  heap = 0,
                                  M = 1,
                                  nitt = 3000,
                                  thin = 10,
                                  burnin = 1000)
          }
        }


        if(tmp_type == "roundedcont"){

          #yet, we don't have a multilevel imputation for rounded incomes
            imp <- imp_roundedcont(y_imp_multi = my_data[, tmp_variable],
                                   X_imp_multi = tmp_X,
                                   intercept_varname = fe$intercept_varname,
                                   M = 1)

        }


        if(tmp_type == "categorical"){
          if(is.null(fe$clID_varname) || !(fe$clID_varname %in% names(my_data)) || tmp_variable != fe$target_varname){

            imp <- imp_cat_cart(y_imp_multi = my_data[, tmp_variable],
                                X_imp_multi = tmp_X)

          }else{

            imp <- imp_cat_multi(y_imp_multi = as.factor(my_data[, tmp_variable]),
                                 X_imp_multi = tmp_X,
                                 Z_imp_multi = tmp_Z,
                                 clID = my_data[, fe$clID_varname],
                                 model_formula = model_formula,
                                 M = 1,
                                 nitt = 3000,
                                 thin = 10,
                                 burnin = 1000)
          }
        }

        if(tmp_type == "ordered_categorical"){
          if(is.null(fe$clID_varname) || !(fe$clID_varname %in% names(my_data)) || tmp_variable != fe$target_varname){
            imp <- imp_orderedcat(data.imp.ordered.categorical = my_data, y.variable.name = tmp_variable,
                                  M = 1)
          }else{
            imp <- imp_orderedcat_multi(y_imp_multi = as.factor(my_data[, tmp_variable]),
                                 X_imp_multi = tmp_X,
                                 Z_imp_multi = tmp_Z,
                                 clID = my_data[, fe$clID_varname],
                                 model_formula = model_formula,
                                 M = 1,
                                 nitt = 3000,
                                 thin = 10,
                                 burnin = 1000)
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

          to_evaluate_numerical <- as.numeric(as.factor(my_data[, tmp_variable]))

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



