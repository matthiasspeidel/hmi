
#help function to do a sample imputation
#' Sample imputation.
#'
#' Function to sample values in a variable from other (observed) values in this variable.
#' So this imputation doesn't use further covariates.
#' @param variable A vector of size \code{n} with missing values.
#' @return A vector of size \code{n} without missing values.
#' @examples
#' sample_imp(c(1, NA, 3, NA, 5))
sample_imp <- function(variable){
  #!!!Es waere sinnvoll hier schon Nebenbedingungen zuzulassen!!!
  ret <- variable
  ret[is.na(ret)] <- sample(size = sum(is.na(variable)),
                             variable[!is.na(variable)], replace = TRUE)
  return(ret)
}


#TODO: Fallunterscheidungen um mit Vectoren, Matrizen oder data.frames umzugehen!!!
#Variable sollte ein Vektor oder factor sein

#' Get the type of variables.
#'
#' Function checks wether a variable is: ...
#' \itemize{
#'  \item continuous (just many numeric values),
#'  \item semicontinuous (many numeric values but more than 5% of the data are 0),
#'  \item a intercept (just the same value for all observations),
#'  \item binary (just two different values - like 0s and 1s or "m" and "f"),
#'  \item categorical (the variable is a factor or has more than 3 different values)
#'  \item orderd categorical (the categorical variable is ordered.)
#'}
#'
#' @param variable A variable from your dataset.
#' @return A character denoting the type of \code{variable}.
get_type <- function(variable){

  if(length(table(variable)) == 1){
    type <- "intercept"
    return(type)
  }

  if(length(table(variable)) == 2){
    type <- "binary"
    #!!!ggfs muessen binaere Variablen noch in 0 und 1 umgewandelt werden.
    #Bsp. wenn die Geschlechtervariable ist bis dahin "m" und "w" ist.
    return(type)
  }

  if(is.vector(variable) || is.factor(variable) || is.array(variable)){
    if(is.numeric(variable)){
      type <- "cont"

      #if the variable is numeric and more than 5% of the variable values are 0,
      #I count this variable as semi-continious.
      #???FROM WHICH RATE ON, ONE CAN SAY, THAT THE VARIABLE IS SEMI-CONT AND NOT
      #ONLY CONT? 1%, 10%, 50%???
      #!!!problematische Zahl frei w?hlbar!!!
      if(sum(variable == 0, na.rm = TRUE)/
         length(variable[!is.na(variable)]) > 0.05){
        type <- "semicont"
      }
      return(type)

    }

    #checking the length of the table to be larger than two is not sufficient
    #to get a real categoriacal variable, because for a (real) continious variable
    #the length of the table is equal to the length of the variable
    #(what will be in general larger than two...)
    #Therefore it will be checked if the table has more than two entries AND
    #wheter the varriable is a factor.
    if(is.factor(variable) || length(table(variable)) > 2){
      type <- "categorical"
      if(is.ordered(variable)){
        type <- "ordered_categorical"
      }

      return(type)
    }

  }else{
    #MS: Wenn es kein Vector oder kein Factor ist, gehe ich davon aus,
    #dass es eine Matrix oder ein data.frame ist.
    #print("R ist komisch!")
    #MS: d.h. ich schliesse aus, dass ein eine Liste mit Datensaetzen etc. ist.
    #MS: Fallunterscheidungen um mit Vectoren, Matrizen oder data.frames umzugehen:
    if(ncol(variable) == 1){
      ret <- get_type(variable[, 1])
      return(ret)
    }
    if(ncol(variable) > 1){
      if(is.matrix(variable)){
        variable <- as.data.frame(variable)
      }
      ret <- unlist(lapply(variable, get_type))
      #MS:?!? Fuer data.frames brauche ich lapply, fuer matrizen apply.
      return(ret)
    }
  }
}



#' Averages the results of the imputation function \code{wrapper}.
#'
#' This function applies the analysis the user is interested in, on all different imputed dataset.
#' Then the results are pooled by simply averaging the results. So the user has to make sure that
#' his analysis produces results with a meaningful average. And furthermore has to accept that no
#' variance is calculated for these parameters.
#' @param imputed_data A \code{mids} object from the \code{wrapper} imputation function.
#' @param analysis A user generated function that he is interested in. See examples.
#' @return A vector with all averaged results.
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
hmi_pool <- function(data, analysis){

  if (!is.mids(data)){
    stop("The data must have class mids.")
  }

  results <-list()

  for (i in 1:data$m) {
    results[[i]] <- analysis(complete(data, i))
  }

  tmp <- simplify2array(results)
  mode(tmp) <- "numeric"
  return(rowMeans(tmp))
}
