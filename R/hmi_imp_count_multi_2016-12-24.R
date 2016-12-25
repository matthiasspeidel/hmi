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



#' The function for hierarchical imputation of contious variables.
#'
#' The function is called by the wrapper.
#' @param y_imp_multi A Vector with the variable to impute.
#' @param X_imp_multi A data.frame with the fixed effects variables.
#' @param Z_imp_multi A data.frame with the random effects variables.
#' @param clID A vector with the cluster ID.
#' @param n.iter An integer defining the number of
#' iterations that should be run in each bunch of iterations.
#' @param M An integer defining the number of imputations that should be made.
#' @param nitt An integer defining number of MCMC iterations (see MCMCglmm).
#' @param thin An integer defining the thinning interval (see MCMCglmm).
#' @param burnin An integer defining the percentage of draws from the gibbs sampler
#' that should be discarded as burn in (see MCMCglmm).
#' @return A n x M matrix. Each column is one of M imputed y-variables.
imp_count_multi <- function(y_imp_multi,
                      X_imp_multi,
                      Z_imp_multi,
                      clID,
                      M = 10,
                      nitt = 3000,
                      thin = 10,
                      burnin = 1000){

  just_cont <- imp_cont_multi(y_imp_multi = y_imp_multi,
                              X_imp_multi = X_imp_multi,
                              Z_imp_multi = Z_imp_multi,
                              clID = clID,
                              M = M,
                              nitt = nitt,
                              thin = thin,
                              burnin =  burnin)

  y_imp <- apply(just_cont, 2, round)
  return(y_imp)
}


# Generate documentation with devtools::document()
# Build package with devtools::build() and devtools::build(binary = TRUE) for zips
