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



#' The main function called by the user. LATER THE USER WILL USE THE WRAPPER.
#'
#' The user has to passes to the function his data.
#' Optionally he pass his analysis model formula so that hmi runs the imputation model
#' in line with his analysis model formula.
#' And of course he can specify some parameters for the imputation routine
#' (like the number of imputations and iterations and the burn in percentage.)
#' The standard usage should be that the user gives his complete dataset
#' and his analysis model. But he also could just give y, X and Z and the cluser ID.
#' @param y_imp_multi A Vector with the variable to impute.
#' @param X_imp_multi A data.frame with the fixed effects variables.
#' @param Z_imp_multi A data.frame with the random effects variables.
#' @param y_variable_name A character naming the variable to impute.
#' @param clID A vector with the cluster ID.
#' @param n.iter An integer defining the number of
#' iterations that should be run in each bunch of iterations.
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
imp_multi <- function(y_imp_multi,
                      X_imp_multi,
                      Z_imp_multi,
                      clID,
                      M = 10,
                      nitt = 3000,
                      thin = 10,
                      burnin = 1000,
                      allowed_max_value = Inf,
                      allowed_max_variable = NULL,
                      allowed_min_value = -Inf,
                      allowed_min_variable = NULL){

  a <- Sys.time()

  # -----------------------------preparing the data ------------------
  # -- standardise the covariates in X (which are numeric and no intercept)
  need_stand <- apply(X_imp_multi, 2, get_type) == "cont"
  X_imp_multi_stand <- X_imp_multi
  X_imp_multi_stand[, need_stand] <- scale(X_imp_multi[, need_stand])
  #X_imp_multi %>% mutate_each_(funs(scale), vars = names(need_stand)[need_stand])
  #generate model.matrix (from the class matrix)
  n <- nrow(X_imp_multi_stand)
  X_model_matrix <- model.matrix(rnorm(n) ~ 0 + ., data = X_imp_multi_stand)
  # Remove ` from the variable names
  colnames(X_model_matrix) <- gsub("`", "", colnames(X_model_matrix))

  # -- standardise the covariates in Z (which are numeric and no intercept)
  need_stand <- apply(Z_imp_multi, 2, get_type) == "cont"
  Z_imp_multi_stand <- Z_imp_multi
  Z_imp_multi_stand[, need_stand] <- scale(Z_imp_multi[, need_stand])

  # MS: If the user wants a fixed intercept, the wrapper function assures that X_imp_multi
  # includes such a variable

  # Get the number of random effects variables
  n.par.rand <- ncol(Z_imp_multi_stand)
  length.alpha <- length(table(clID)) * n.par.rand


  # -------------- calling the gibbs sampler to get imputation parameters----

  tmp <- data.frame(target = y_imp_multi)
  xnames <- paste("X", 1:ncol(X_model_matrix), sep = "")
  znames <- paste("Z", 1:ncol(Z_imp_multi_stand), sep = "")
  tmp[, xnames] <- X_model_matrix
  tmp[, znames] <- Z_imp_multi_stand
  tmp[, "ClID"] <- clID

  fixformula <- formula(paste("target~", paste(xnames, collapse = "+"), "- 1", sep = ""))
  randformula <- as.formula(paste("~us(", paste(znames, collapse = "+"), "):ClID", sep = ""))


  prior <- list(R = list(V = 1e-07, nu = -2),
                G = list(G1 = list(V = diag(2), nu = 2))) #PRIORIS NOCH ANPASSEN
  #VOR ALLEM AUCH FLEXIBEL MACHEN diag(ncol(Z_imp_multi_stand)) etc.
  MCMCglmm_draws <- MCMCglmm::MCMCglmm(fixformula, random = randformula, data = tmp,
                             verbose = FALSE, pr = TRUE, prior = prior,
                             saveX = TRUE, saveZ = TRUE,
                             nitt = 3000,
                             thin = 10,
                             burnin = 1000)

  pointdraws <- MCMCglmm_draws$Sol
  xdraws <- pointdraws[, 1:ncol(X_model_matrix), drop = FALSE]
  zdraws <- pointdraws[, ncol(X_model_matrix) + 1:length.alpha, drop = FALSE]
  variancedraws <- MCMCglmm_draws$VCV
  #! Die letzte Spalte beinhaltet die VARIANZ (nicht die Standardabweichung)
  # Der Residuen!
  number_of_draws <- nrow(pointdraws)
  select.record <- sample(1:number_of_draws, M, replace = TRUE)

  # -------------------- drawing samples with the parameters from the gibbs sampler --------
  y.imp <- array(NA, dim = c(n, M))
  ###start imputation
  for (j in 1:M){

    rand.eff.imp <- matrix(zdraws[select.record[j],],
                           ncol = n.par.rand)

    fix.eff.imp <- matrix(xdraws[select.record[j], ], nrow = ncol(X_model_matrix))

    sigma.y.imp <- sqrt(variancedraws[select.record[j], ncol(variancedraws)])

    y.temp <- rnorm(n, X_model_matrix %*% fix.eff.imp +
                      apply(Z_imp_multi_stand * rand.eff.imp[clID,], 1, sum), sigma.y.imp)

    y.imp[, j] <- ifelse(is.na(y_imp_multi), y.temp, y_imp_multi)
  }
  b <- Sys.time()
  print(b - a)
  # --------- returning the imputed data --------------
  return(y.imp)
}


# Generate documentation with devtools::document()
# Build package with devtools::build() and devtools::build(binary = TRUE) for zips
