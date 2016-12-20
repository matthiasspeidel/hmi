# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'
#
# MS: All other shortcuts are shown when pressing 'Alt + Shift + K'
# MS: Hinweis: es wird zu jeder R-Datei im Projekordner\R, die eine Documentation hat,
#Mit devtools::document() auch eine .RD Datei erstellt.


#library(MASS) 'MS: DURCH BESSEREN AUFRUF require (etc) ersetzen
#library(coda) #is needed for gelman.diag()
#library(stats)

#########three different parts:

#########Gibbs sampler functions: contains all functions to draw from
#########				    the conditional distributions
#########Gibbs sampler: 	    the actual Gibbs sampler function
#########Imputation function:     the function that prepares the data before
#########				   calling the Gibbs sampler and uses
#########			          the parameters from the Gibbs for imputation







#################################################
#########Imputation function#####################
#################################################



#' The function for hierarchical imputation of categorical variables.
#'
#' The function is called by the wrapper.
#' @param y_imp_multi A Vector with the variable to impute.
#' @param X_imp_multi A data.frame with the fixed effects variables.
#' @param Z_imp_multi A data.frame with the random effects variables.
#' @param clID A vector with the cluster ID.
#' @param model_formula A \code{\link[stats]{formula}} used for the analysis model.
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
imp_cat_multi <- function(y_imp_multi,
                      X_imp_multi,
                      Z_imp_multi,
                      clID,
                      model_formula,
                      M = 10,
                      nitt = 3000,
                      thin = 10,
                      burnin = 1000,
                      allowed_max_value = Inf,
                      allowed_max_variable = NULL,
                      allowed_min_value = -Inf,
                      allowed_min_variable = NULL){

  if(!is.factor(y_imp_multi)){
    warning("We suggest to make all your categorical variables to be a factor.")
    y_imp_multi <- as.factor(y_imp_multi)
  }

  # -----------------------------preparing the data ------------------
  # -- standardise the covariates in X (which are numeric and no intercept)
  need_stand <- apply(X_imp_multi, 2, get_type) == "cont"
  X_imp_multi_stand <- X_imp_multi
  #X_imp_multi_stand[, need_stand] <- scale(X_imp_multi[, need_stand])
  #X_imp_multi %>% mutate_each_(funs(scale), vars = names(need_stand)[need_stand])
  #generate model.matrix (from the class matrix)
  n <- nrow(X_imp_multi_stand)
  tmp0 <- data.frame(y_tmp = rnorm(n), X_imp_multi_stand, cl.id = clID)
  tmp0_formula <- formula(paste("y_tmp ~ 0 + ", paste(colnames(X_imp_multi_stand), collapse = "+")))

  X_model_matrix <- model.matrix(tmp0_formula, data = tmp0)
  #!!! BESSER LOESEN!!!
    #alt: model.matrix(rnorm(n) ~ 0 + ., data = X_imp_multi_stand)
  # Remove ` from the variable names
  colnames(X_model_matrix) <- gsub("`", "", colnames(X_model_matrix))



  # -- standardise the covariates in Z (which are numeric and no intercept)
  need_stand <- apply(Z_imp_multi, 2, get_type) == "cont"
  Z_imp_multi_stand <- Z_imp_multi
  #Z_imp_multi_stand[, need_stand] <- scale(Z_imp_multi[, need_stand])

  # MS: If the user wants a fixed intercept, the wrapper function assures that X_imp_multi
  # includes such a variable


  number_clusters <- length(table(clID))


  # -------------- calling the gibbs sampler to get imputation parameters----

  tmp <- data.frame(target = y_imp_multi)
  xnames <- paste("X", 1:ncol(X_model_matrix), sep = "")
  znames <- paste("Z", 1:ncol(Z_imp_multi_stand), sep = "")
  tmp[, xnames] <- X_model_matrix[, 1:ncol(X_model_matrix)]
  # note: the [] part is only relevant if ncol(X_model_matrix) == 1
  # (see http://stackoverflow.com/questions/40872034/what-happens-when-a-data-frame-gets-new-columns)
  tmp[, znames] <- Z_imp_multi_stand[, 1:ncol(Z_imp_multi_stand)]
  tmp[, "ClID"] <- clID

  fixformula <- formula(paste("target~ - 1 + ",
                                paste(xnames, ":trait", sep = "", collapse = "+"), sep = ""))

  randformula <- formula(paste("~us( - 1 + ", paste(znames, ":trait", sep = "", collapse = "+"),
                                 "):ClID", sep = ""))


  number_fix_parameters <- ncol(X_model_matrix) * (J-1)

  # Get the number of random effects variables
  number_random_effects <- length(znames)
  number_random_parameters <- number_random_effects * (J - 1)
  #Fix residual variance R at 1
  # cf. http://stats.stackexchange.com/questions/32994/what-are-r-structure-g-structure-in-a-glmm
  #prior <- list(R = list(V = diag(2), nu = 0.002, fix = 1),
    #            G = list(G1 = list(V = diag(2), nu = -1,
      #                             alpha.mu = c(0,0), alpha.V = diag(2) *25^2))) #PRIORIS NOCH ANPASSEN#G1 = list(V = diag(n.par.rand), nu = 0.002)
  J <- length(table(y_imp_multi)) #number of categories
  J_matrix <- array(1, dim = c(J, J) -1) # matrix of ones
  I_matrix <- diag(J - 1) #identiy matrix

  IJ<- (I_matrix + J_matrix)/J # see Hadfields Course notes p 97
  prior <- list(R = list(V = IJ, fix = TRUE),
                G = list(G1 = list(V = diag(number_random_parameters), nu = 2)))
  # Wenn ich es richtig verstehe, ginge auch V = diag(2) und n = -1 oder n = -2
  # fuer flache priors auf ]0, inf] bzw. nicht informative priors.
  # (mir ist also nicht klar wo der unterschied ist).
  # ist nicht jede flache prio auch non-informative?
  MCMCglmm_draws <- MCMCglmm::MCMCglmm(fixed = fixformula,
                                       random = randformula,#~us(trait):ClID,
                                       rcov = ~us(trait):units,
                                       data = tmp,
                                       family = "categorical",
                                       verbose = FALSE, pr = TRUE,
                                       prior = prior,
                                       saveX = TRUE, saveZ = TRUE,
                                       nitt = 1e5,
                                       thin = 1e3,
                                       burnin = 1e4)

  # ??? Maybe a correction helps ???
  # see http://stats.stackexchange.com/questions/32994/what-are-r-structure-g-structure-in-a-glmm
 # k <- ((16*sqrt(3))/(15*pi))^2

  pointdraws <- MCMCglmm_draws$Sol #/sqrt(1 + k)
  #hier müssen jetzt für jede kategorie die regressionsparameter extrahiert werden.

  xdraws <- pointdraws[, 1:number_fix_parameters, drop = FALSE]
  zdraws <- pointdraws[, number_fix_parameters +
                         1:(number_random_parameters * number_clusters), drop = FALSE]

  #! Die letzte Spalte beinhaltet die VARIANZ (nicht die Standardabweichung)
  # Der Residuen!
  number_of_draws <- nrow(pointdraws)
  select.record <- sample(1:number_of_draws, M, replace = TRUE)

  # -------------------- drawing samples with the parameters from the gibbs sampler --------
#################################### BEGINN new stuff ################
  #now generate new P(Y = A|x * beta) = x*beta/(1+ sum(exp(x*beta))) etc.

  #set up random intercepts and slopes


  exp_beta <- array(dim = c(n, J - 1))

  fix_eff_sample <- xdraws[select.record[1], ]
  rand_eff_sample <- zdraws[select.record[1], ]
  for(k in 1:ncol(exp_beta)){

    rand_betas <- matrix(rand_eff_sample[(k - 1) * number_clusters +
                     c(1:number_clusters, 1:number_clusters  + (J - 1) * number_clusters)], nrow = number_clusters)

    exp_beta[, k] <-
      as.matrix(exp(as.matrix(tmp[, xnames, drop = FALSE]) %*% fix_eff_sample[k + c(0, J - 1)] +#fix effects
      rowSums(tmp[, znames, drop = FALSE] * rand_betas[clID, , drop = FALSE])))#random effects

  }


  ####################### end new stuff #############


  y_imp <- data.frame(matrix(nrow = n, ncol = M))
  ###start imputation
  for (j in 1:M){

    y_temp <- array(dim = n)

    #VERBESSERUNGSPOTENTIAL: VEKTOWERTIG MACHEN
    for(i in 1:n){

      #!!!MAKE SHURE THAT THE REFERENCE CATEGORY IS ON THE CORRECT POSITION
      mytmp_prob <- exp_beta[i, ]/(1 + sum(exp_beta[i, ]))
      my_prob <- c(1 - sum(mytmp_prob), mytmp_prob) #the first category is the reference category in MCMCglmm

      #c(mytmp_prob, 1 - sum(mytmp_prob))

      y_temp[i] <- sample(levels(y_imp_multi), size = 1, replace = TRUE, prob = my_prob)

    }

    y_imp[, j] <- as.factor(ifelse(is.na(y_imp_multi), y_temp, as.character(y_imp_multi)))
  }
  colnames(y_imp) <- paste("imp", 1:M, sep = "_")

  # --------- returning the imputed data --------------
  return(y_imp)
}


# Generate documentation with devtools::document()
# Build package with devtools::build() and devtools::build(binary = TRUE) for zips
