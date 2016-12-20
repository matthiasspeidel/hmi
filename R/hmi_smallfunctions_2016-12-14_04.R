
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
#'  \item continuous (numeric values),
#'  \item semicontinuous (numeric values with more than 5% of them are 0),
#'  \item rounded continuous (if continuous values are rounded to the closest multiple of 5, 10, 50, 100, 500 or 1000. We see this to be the case if more than 50 % are devisable by 5)
#'  \item a intercept (the same value for all observations),
#'  \item binary (two different values - like 0s and 1s or "m" and "f"),
#'  \item categorical (the variable is a factor or has more than 3 different values)
#'  \item orderd categorical (the categorical variable is ordered.)
#'}
#'
#' @param variable A variable (vector) from your dataset.
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

      if(sum(variable %% 5 == 0, na.rm = TRUE)/
        length(variable[!is.na(variable)]) > 0.5){
        type <- "roundedcont"
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


#' calculate the likelihood contribution of the data
#'
#' Needed for roundedcont
#' @param para is the vector of paramerers and according to them, the value negloglik2 returns shall be
#' maximized. So the task is: find the setting of para maximizing negloglik2(para, ...)
#' The starting values are c(kstart, betastart2, gammastart, sigmastart2)
#' (the 6 treshholds (or "cutting points") for the latent variable behind the rounding degree,
#' the regression parameters explaining the logged income)
#' @param X the data.frame of covariates
#' @param y_in_negloglik the target variable
#' @param myp This vector is the indicator of the (highes possible) rounding degree for an observation.
#' This parameter comes directly from the data.
#' @param mean.log.inc Is the integer (!!!WHEN I WRITE "INTEGER" I DONT MEAN ONLY -2, -1, 0, 1, 2,
#' ... BUT ANY SINGLE NUMBER LIKE -1.432423, pi, 2.3444, 1, ... I USE THIS NOTATION AS "NUMBER" OFTEN MEANS
#' THE RESULT WHEN COUNTING SOMETHING: THE NUMBERS OF SIMULATION RUNS OR THE NUMBER OF OBSERVATIONS, SO ONLY INTEGERS.)
#' with the value of the mean of the logarithm of the target value
#' @param sd.log.inc Is the integer with the value equal to the standard deviation of the logarithm of the target variable
#' @return An integer equal to the (sum of the) negative log-likelihood contributions (of the observations)
negloglik2 <- function(para, X, y_in_negloglik, myp, mean.log.inc, sd.log.inc){

  #MS: jedes k ist ein unbekannter "treshold value" (Paper Seite 12)
  #MS: und jenachdem zwischen welchen k die Variable G (die bedingt auf die anderen Kovariablen und Y normalverteilt ist) liegt
  #MS: liegt ein anderer Rundungsprozess vor.
  k0 <- para[1]
  k1 <- para[2]
  k2 <- para[3]
  k3 <- para[4]
  k4 <- para[5]
  k5 <- para[6]

  #the regression coefficients beta defining mean2 (the expected value of log(Y), see eq(2) in Drechsler and Kiesl, 2014)
  #also the appear in mean1 (the expected value of G)
  reg.coefs <- matrix(para[7:(length(para) - 2)], ncol = 1)

  gamma1 <- para[length(para) - 1] #the gamma1 from eq (3)

  sigmax <- para[length(para)] # the variance for log(Y) (see Drechsler, Kiesel, 2014) equation (3) page 12

  if (k1 < k0) return(1e50)        # make sure that threshholds
  if (k2 < k1) return(1e50)        # are in increasing order
  if (k3 < k2) return(1e50)        #
  if (k4 < k3) return(1e50)        #
  if (k5 < k4) return(1e50)        #
  #MS: kuerzer  waere: if(is.unsorted(par[1:6], strictly = TRUE)) return(1e50)
  #MS: Die funktion prueft ob par[1:6] (streng) monoton steigend ist.

  n <- nrow(X)
  result <- rep(NA, n)

  # mean and covariance of the joint normal of x and g ,
  # following Rubin/Heitjan

  #the mean for G (see Drechsler, Kiesel, 2014) equation (2) page 12
  mean1 <- gamma1 * (as.matrix(X) %*% reg.coefs)

  #?By not using Z, we assume, that all variables that might be in Z are already in X?

  #the mean for log(Y) (see Drechsler, Kiesel, 2014) equation (2) page 12
  mean2 <- as.matrix(X) %*% reg.coefs #MS: ist also grob gesprochen y_dach

  #the covariance matrix for log(Y) and G (see Drechsler, Kiesel, 2014) equation (3) page 12
  sigma <- matrix(c(1 + gamma1^2 * sigmax^2,
                    gamma1 * sigmax^2,
                    gamma1 * sigmax^2,
                    sigmax^2), nrow = 2)


  #MS: allerdings ist tau0 auf 1 festgesetzt (s. 3.2.2) und hier ist G,Y angegebene,
  #MS: im Paper aber Y,G.
 # print(cov2cor(sigma)[1, 2])
  rho <- max(min(cov2cor(sigma)[1, 2], 1), -1)  #MS: kuerzer als sigma[1, 2]/sqrt(sigma[1, 1] * sigma[2, 2])
#print(rho) #wird manchmal schlagartig 1
  #print(paste("rho after correction:", rho))
  sigma1 <- sqrt(sigma[1, 1])
  sigma2 <- sqrt(sigma[2, 2])

  #MS: Extrahieren der beobachteteten Einkommen (wird spaeter fuer pbivnormX() gebraucht)
  missind_negloglik <- is.na(y_in_negloglik)
  inc.obs <- y_in_negloglik[!missind_negloglik]

  # a1 are a7 the "components" of the log-liklihood contributions
  #
  # If the income is not divisable by 5, only a1 is relevant
  # if the income is divisable by 5, but not by 10, a1 + a2 are relevant
  # etc.

  #get the value of the standard bivariate normal cumulative density funtion
  #the better k0, k1, ..., sigma and rho are found, the better a1, a2, ...

  a1 <- pbivnormX(x =  c(k0 - mean1[!missind_negloglik])/sigma1,
                  y =  ((log(inc.obs + 0.5) - mean.log.inc)/sd.log.inc - mean2[!missind_negloglik])/sigma2,
                  rho = rho
  ) -
    pbivnormX(x =  c(k0 - mean1[!missind_negloglik])/sigma1,
              y =  ((log(pmax(inc.obs - 0.5, 0)) - mean.log.inc)/sd.log.inc - mean2[!missind_negloglik])/sigma2,
              rho = rho
    )

  a1 <- pmax(1e-100, abs(a1)) # reason: if a1 is close to 0, a1 can be 0
  # due to computational imprecision

  a2 <- components(ki = k0, kj = k1,
                   mean1_obs = mean1[!missind_negloglik],
                   mean2_obs = mean2[!missind_negloglik],
                   sigma1 = sigma1,
                   sigma2 = sigma2,
                   rho = rho,
                   inc_obs = inc.obs,
                   mean.log.inc = mean.log.inc,
                   sd.log.inc = sd.log.inc,
                   half.divisor = 2.5)

  a3 <- components(ki = k1, kj = k2,
                   mean1_obs = mean1[!missind_negloglik],
                   mean2_obs = mean2[!missind_negloglik],
                   sigma1 = sigma1,
                   sigma2 = sigma2,
                   rho = rho,
                   inc_obs = inc.obs,
                   mean.log.inc = mean.log.inc,
                   sd.log.inc = sd.log.inc,
                   half.divisor = 5)

  a4 <- components(ki = k2, kj = k3,
                   mean1_obs = mean1[!missind_negloglik],
                   mean2_obs = mean2[!missind_negloglik],
                   sigma1 = sigma1,
                   sigma2 = sigma2,
                   rho = rho,
                   inc_obs = inc.obs,
                   mean.log.inc = mean.log.inc,
                   sd.log.inc = sd.log.inc,
                   half.divisor = 25)

  a5 <- components(ki = k3, kj = k4,
                   mean1_obs = mean1[!missind_negloglik],
                   mean2_obs = mean2[!missind_negloglik],
                   sigma1 = sigma1,
                   sigma2 = sigma2,
                   rho = rho,
                   inc_obs = inc.obs,
                   mean.log.inc = mean.log.inc,
                   sd.log.inc = sd.log.inc,
                   half.divisor = 50)

  a6 <- components(ki = k4, kj = k5,
                   mean1_obs = mean1[!missind_negloglik],
                   mean2_obs = mean2[!missind_negloglik],
                   sigma1 = sigma1,
                   sigma2 = sigma2,
                   rho = rho,
                   inc_obs = inc.obs,
                   mean.log.inc = mean.log.inc,
                   sd.log.inc = sd.log.inc,
                   half.divisor = 250)


  a7 <- pnorm((log(inc.obs + 500) - mean.log.inc)/sd.log.inc, mean2[!missind_negloglik], sigma2) -

    pbivnormX( x =  c(k5 - mean1[!missind_negloglik])/sigma1,
               y =  ((log(inc.obs + 500) - mean.log.inc)/sd.log.inc - mean2[!missind_negloglik])/sigma2,
               rho = rho
    ) -

    pnorm((log(pmax(inc.obs - 500, 0)) - mean.log.inc)/sd.log.inc, mean2[!missind_negloglik], sigma2)  +

    pbivnormX( x =  c(k5 - mean1[!missind_negloglik])/sigma1,
               y =  ((log(pmax(inc.obs - 500, 0)) - mean.log.inc)/sd.log.inc - mean2[!missind_negloglik])/sigma2,
               rho = rho
    )

  a2 <- pmax(0, a2)#MS: replacing negative values by 0.
  a3 <- pmax(0, a3)
  a4 <- pmax(0, a4)
  a5 <- pmax(0, a5)
  a6 <- pmax(0, a6)
  a7 <- pmax(0, a7)

  #MS: Q: Was ist myp?
  #MS: A: Der Indikator, in welchem Intervall man sich befindet.
  result_try <- cbind(a1 , a2 * (myp[!missind_negloglik] >= 1) , a3 * (myp[!missind_negloglik] >= 2),
                      a4 * (myp[!missind_negloglik] >= 3) , a5 * (myp[!missind_negloglik] >= 4),
                      a6 * (myp[!missind_negloglik] >= 5), a7 * (myp[!missind_negloglik] == 6))
  #MS: result.try ist eine 8179 x 7 Matrix in der es fuer jede Beobachtung mit praeziser Einkommensangabe
  #MS: eine eigene Zeile gibt. In jeder dieser Zeile beeinhaltet die erste Spalte den Wert von a1,
  #MS: In der zweiten Spalte steht der Wert von a2, sofern diese Person ihr Einkommen mindestens zum
  #MS: zweiten Rundungsgrad (hier ohne Rest teilbar durch 5) angegeben hat. Ansonsten ist der Wert 0.

  #es ist also der logarithmierte beitrag zur likelihood (oder anders ausgedr?ckt:
  # der beitrag zur log-likelihood) f?r die Beobachtungen.
  #F?r alle Grade zu denen gerundet werden k?nnte, wird der Likelihood beitrag ausgegeben.
  #Und diese Beitr?ge werden aufsummiert (eq (5) in Drechsler, Kiesl, Speidel, 2015)
  result <- log(rowSums(result_try, na.rm = TRUE))

  ret <- sum(result)
  return(-ret) #notice the minus
}

###########################################################################
## corrected version of pbivnorm from the same package#####################
###########################################################################
#MS: Function pbivnorm() from package "pbivnorm" with small modifications.
#MS: It "[calculates] probabilities from the CDF [Cumulative distribution function]
#MS: of a standard bivariate normal distribution".
#MS: x ist dabei die obere Integrationsgrenze fuer die erste Variable und y fuer die zweite
#MS: (also das x und y in F(x, y) = P(X <= x & Y <= y)) (?).
#MS: Q: wozu braucht man das rho, wenn es doch standard bivarate normalverteilt ist?
#MS: A: Offenbar sagt das "standard" in "Standard bivariate normal" nur, dass die varianzen 1 sind.
#MS: Die Korrelation ist beliebig.

#' calculate probabilities from the cumulative distribution function of a standard bivariate normal distribution
#'
#' A modfied version of pbivnorm() from package "pbivnorm"
#' @param x the value of the first covariate
#' @param y the value of the second covariate
#' @param rho the correlation between the two covariates.
pbivnormX <- function (x, y, rho = 0) {
  if  (is.matrix(x)) {
    #MS: if(2) print("geht in die Schleife")
    #MS: was passiert hier? Die Laenge der Dimension von x kann ja nicht TRUE sein
    #MS: sie ist hoechstens 1, ich wuesste aber nicht in welchem Setting bzw. 0 wenn x keine Matrix sondern eine Zahl ist.
    #MS: Vorschlag: if(!is.matrix(x))
    if (ncol(x) != 2) #MS: also wenn es keine Matrix ist... wie kann ich dann die Anzahl der Zeilen checken?
      stop("'x' must have two columns if specified as a matrix")
    if (!missing(y) && !is.null(y)) #
      warning("'x' was specified as a matrix, so 'y' will be ignored")
    y <- x[, 2]
    x <- x[, 1]
  }
  if (any(abs(rho) > 1)) stop("'rho' must be a valid correlation (-1 <= rho <= 1)")

  if (length(x) != length(y)) stop("'x' and 'y' must have same length")

  x <- replace(x, x == Inf, 1e200)      #  function adjusted here
  x <- replace(x, x == -Inf, -1e200)    #  function adjusted here
  y <- replace(y, y == Inf, 1e200)      #  function adjusted here
  y <- replace(y, y == -Inf, -1e200)    #  function adjusted here

  correl <- numeric(length(x))
  correl[] <- rho
  #MS: Q: Warum nicht gleich correl <- rho?
  #MS: A: Weil sonst correl nicht mehr bspw. 10 Eintraege mit dem Wert von rho haette, sondern nur noch ein Skalar waere
  lower <- as.double(c(0, 0))
  infin <- as.integer(c(0, 0))
  uppera <- as.double(x)
  upperb <- as.double(y)
  lt <- as.integer(length(x))
  prob <- double(lt)
  correl <- as.double(correl)
  ans <- .Fortran("PBIVNORM", prob, lower, uppera, upperb,
                  infin, correl, lt, PACKAGE = "pbivnorm")[[1]]
  return(ans)
}

#' Function to get the likelihood contribution of different roudinding degrees
#'
#' See equation (5) in Drechsler, Kiesel and Speidel (2015)
#' @param ki An integer for the i-th treshold (or "cutpoint")
#' @param kj An integer for the j-th treshold (or "cutpoint") (ki < kj)
#' @param mean1_obs A vector with the expected value of G for the observed data
#' @param mean2_obs A vector with the expected value of log(Y)for the observed data
#' @param sigma1 An integer for the variance of the G
#' @param rho An Integer from [-1, 1] specifying the correlation between G and log(Y)
#' @param inc_obs The vector of the target variable (with all NAs removed)
#' @param mean.log.inc An integer specifying the mean of the logarithm of the target variable
#' @param sd.log.inc An integer specifying the standard deviation of the logarithm of the target variable
#' @param half.divisor An integer needed to find the bounds of possible rounding.
#' If rounding happens to the closest multiple of 500, the half.divisor is 250.
#' @return A vector with the contribution to the likelihood of every individual with an observed target variable value
components <- function(ki, kj, mean1_obs, mean2_obs, sigma1, sigma2, rho, inc_obs,
                       mean.log.inc, sd.log.inc, half.divisor){

  upper_inner <- ((log(inc_obs + half.divisor) - mean.log.inc)/sd.log.inc - mean2_obs)/sigma2
  lower_inner <- ((log(pmax(inc_obs - half.divisor, 0)) - mean.log.inc)/sd.log.inc - mean2_obs)/sigma2
  lower_outer <- c(ki - mean1_obs)/sigma1
  upper_outer <- c(kj - mean1_obs)/sigma1
  ret <- doubleintegral(lower_inner = lower_inner, upper_inner = upper_inner,
                        lower_outer = lower_outer, upper_outer = upper_outer,
                        cdf = pbivnormX, rho = rho)

  return(ret)
}

#' Function to calculate double integrals
#'
#' See eq (5) Drechsler, Kiesel and Speidel (2015)
#' @param lower_inner The vector for the lower bound for the inner integral
#' @param upper_inner The vector for the upper bound for the inner integral
#' @param lower_outer The vector for the lower bound for the outer integral
#' @param upper_outer The vector for the upper bound for the outer integral
#' @param cdf the cumulative density function (from the class "function")
#' @return a vector with the value of the double integral for each observation (with an observed target variable)
doubleintegral <- function(lower_inner, upper_inner, lower_outer, upper_outer,
  cdf, ...){
  ret <- (
    cdf(x = upper_outer, y =  upper_inner, ...) -
      cdf(x = upper_outer, y =  lower_inner, ...)
  ) - (
    cdf(x = lower_outer, y =  upper_inner, ...)-
      cdf(x = lower_outer, y =  lower_inner, ...)
  )
  return(ret)
}


#' Function to extract the different elements of a formula
#'
#' The function shall find: the target variable, the variables with the fixed effects,
#' if there is a cluster ID: this and the random effects variables.
#' The names of the fixed and random intercepts variable are explicitly labelled.
#' @param model_formula A formula (from class formula)
#' @param constant_variables A Boolean-vector of length equal to the number of columns in the data
#' specifying whether a variable is a constant variable (eg. an intercept variable)
#' @param variable_names_in_data A character-vector with the column names of the data.
#' @return a list with the names of the target variable, the intercept variable,
#' the fixed and random effects covariates and the cluster id variable.
extract_varnames <- function(model_formula = NULL, constant_variables,
                             variable_names_in_data){


  # Set up default values for key variables like the target variable, the clusterID,
  # and the random covariates
  target_varname <- NULL
  #extract key variables from formula
  if(!is.null(model_formula)){

    # if it was a character string, make it a formula
    if(class(model_formula) == "character"){
      warning("We need your model_formula to be a formula (and not a character).
              So we changed it automatically, but it would be better if you do it.")

      model_formula <- formula(model_formula)
    }

    variable_names_in_formula <- all.vars(terms(model_formula))


    # ----------Target variable
    target_varname_full <- as.character(model_formula)[2]

    #check if the formula contains any functions like "log(y)" in log(y) ~ x.
    target_varname <- gsub(".*\\((.*)\\).*", "\\1", target_varname_full)

    if(target_varname != target_varname_full){
      warning(paste("For the tool to work saver we suggest to do the transformation -->", target_varname_full,
                    "<-- before running the imputation."))
    }

    # -------- Cluster ID
    #clID_varname <- as.character(lme4::findbars(model_formula)[[1]])[3]

    clID_varname <- sapply(lme4::findbars(model_formula), function(x) as.character(x)[3])

    #check if there is a cluster variable specified:...
    if(any(is.na(clID_varname)) | length(clID_varname) == 0){
      #... if it is not, we don't have a random effects model
      clID_varname <- ""
      random_intercept_exists <- FALSE
      randomeffects_varname <- ""

    }else{#if the cluster id is not NA


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

    }

    # ---------X variables

    fixedeffects_varname <- all.vars(delete.response(terms(lme4::nobars(model_formula))))
    # --check if intercept is present
    fixed_intercept_exists <- attributes(delete.response(terms(lme4::nobars(model_formula))))$intercept == 1

    # --remove cluster ID from the X-variables (but it shouldn't be there)
    fixedeffects_varname <- fixedeffects_varname[ ! fixedeffects_varname %in% clID_varname]


    if("0" %in% fixedeffects_varname){
      fixed_intercept_exists <- FALSE
      #exclude fixed effects called "0"
      fixedeffects_varname <- fixedeffects_varname[fixedeffects_varname != "0"]
    }

    ## Handling of a intercept variable.
    # MS: for the imputation functions it is important,
    # that constant variables are handled correctly


    # make shure, the data will have a "Intercept" variable, if the model_formula says to model an intercept.
    intercept_varname <- NULL
    if(fixed_intercept_exists | random_intercept_exists | sum(constant_variables) == 1){

      # give intercept_varname a default value.
      intercept_varname <- "Intercept"
      #MS:?!? If there is already an intercept, consider renaming it to "(Intercept)"
      # for reasons of consitency with further code.
      #MS: If there is no intercept yet in the data, generate one

      if(sum(constant_variables) == 1){
        intercept_varname <- variable_names_in_data[constant_variables]
      }


      if(fixed_intercept_exists){
        # MS: man muss sicher stellen, dass der alte intercept name nicht noch zusätzlich auftaucht
        # bspw war bei den random effects "1" ein Element von randomeffects_varname.
        # der Intercept war aber "Intercept", so dass die Schleife Stand 2016-07-17
        # "Intercept" "1" "X"  ausgab.
        fixedeffects_varname <- c(intercept_varname, fixedeffects_varname)

        # If the intercept_varname is not "1"...
        if(intercept_varname != "1"){
          #... then every fixedeffects_varname that is "1" has to be removed
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


  }else{ # if model_formula is NULL
    # in this case impute every variable based on the others in a single level framework
    clID_varname <- NULL
    fixedeffects_varname <- colnames(my_data)
    randomeffects_varname <- character()
  }

  ret <- list(target_varname = target_varname,
              intercept_varname = intercept_varname,
              fixedeffects_varname = fixedeffects_varname,
              randomeffects_varname = randomeffects_varname,
              clID_varname = clID_varname)
  return(ret)
}
