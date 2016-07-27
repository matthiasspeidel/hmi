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
#' @param n.chains An integer defining the number of Markov chains to be made.
#' @param burn.in A numeric between 0 and 1 defining the percentage of draws from the gibbs sampler
#' that should be discarded as burn in.
#' @param max.iter An integer defining the maximum number of
#' iterations that should be run in total.
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
                      n.iter = 100, M = 10, n.chains = 3, burn.in = 1/3,
                      max.iter = 5000,
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


  #------- beginn testfield
  prior <- list(R = list(V = 1e-07, nu = -2),
                G = list(G1 = list(V = diag(2), nu = 2)))
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
  MCMCglmm_draws <- MCMCglmm(fixformula, random = randformula, data = tmp,
                             verbose = FALSE, pr = TRUE, prior = prior,
                             saveX = TRUE, saveZ = TRUE,
                             nitt = 3000,
                             thin = 10,
                             burnin = 1000)

  pointdraws <- MCMCglmm_draws$Sol
  xdraws <- pointdraws[, 1:ncol(X_model_matrix), drop = FALSE]
  zdraws <- pointdraws[, ncol(X_model_matrix) + 1:length.alpha, drop = FALSE]
  variancedraws <- MCMCglmm_draws$VCV
  number_of_draws <- nrow(pointdraws)
  select.record <- sample(1:number_of_draws, M, replace = TRUE)


  #---------end testfield

 # pars <- gibbs.coef(y_gibbs = y_imp_multi, X_gibbs = X_model_matrix, Z_gibbs = Z_imp_multi_stand,
##                     clID = clID,
 #                    n.iter = n.iter, M = M, n.chains = n.chains, burn.in = burn.in,
#                     max.iter = max.iter)

  ###sample m values for each parameter
 # select.record <- sample(1:dim(pars)[2], M, replace = TRUE)
 # select.chain <- sample(1:dim(pars)[1], M, replace = TRUE)

  # -------------------- drawing samples with the parameters from the gibbs sampler --------
  y.imp <- array(NA, dim = c(n, M))
  ###start imputation
  for (j in 1:M){

    rand.eff.imp <- matrix(zdraws[select.record[j],],
                           ncol = n.par.rand)

    fix.eff.imp <- matrix(xdraws[select.record[j], ], nrow = ncol(X_model_matrix))

    sigma.y.imp <- variancedraws[select.record[j], ncol(variancedraws)]

    y.temp <- rnorm(n, X_model_matrix %*% fix.eff.imp +
                      apply(Z_imp_multi_stand * rand.eff.imp[clID,], 1, sum), sigma.y.imp)

    y.imp[, j] <- ifelse(is.na(y_imp_multi), y.temp, y_imp_multi)
  }
  b <- Sys.time()
  print(b - a)
  # --------- returning the imputed data --------------
  return(y.imp)
}


############################################################
######Gibbs sampler for random coefficient model############
############################################################

#' Gibbs sampler for random coefficient model
#'
#' It generates the Markov chains of the imputation parameters by drawing from their conditional distributions
#' until convergence
#' @param y_gibbs A vector or data.frame with \code{ncol = 1} containing the target variable with the missing values.
#' @param X_gibbs A matrix containing the covariates influencing \code{y} via fixed effects.
#' If rows with missing values in \code{X} should also be imputed, put all your variables in a data.frame (or matrix)
#' @param Z_gibbs A data.frame containing the covariates influencing \code{y} via random effects
#' @param clID A factor (should come as data.frame or vector) containing the cluster IDs.
#' @param n.iter An integer defining the number of
#' iterations that should be run in each bunch of iterations.
#' @param max.iter An integer defining the maximum number of
#' iterations that should be run in total.
#' @param M An integer defining the number of imputations that should be made.
#' @param n.chains An integer defining the number of Markov chains to be made.
#' @param burn.in A numeric between 0 and 1 defining the percentage of draws from the gibbs sampler
#' that should be discarded as burn in.
#' @return It returns a multidimensional vector with the Markov chains
#' containing the imputation parameters needed for \code{imp_multi}.
gibbs.coef <- function(y_gibbs, X_gibbs, Z_gibbs, clID, n.iter = 100, M = 10,
                       n.chains = 3, burn.in = 1/3, max.iter = 5000){

  ##calculate quantities that don't change in the Gibbs sampler
  y_obs <- y_gibbs[!is.na(y_gibbs)] ###only observed part of y
  X_obs <- X_gibbs[!is.na(y_gibbs), , drop = FALSE]
  Z_obs <- Z_gibbs[!is.na(y_gibbs), , drop = FALSE]
  clID_obs <- clID[!is.na(y_gibbs)]
  clID_obs <- droplevels(clID_obs)

  n.obs <- sum(!is.na(y_gibbs))
  nclass <- length(table(clID))
  n.par.rand <- ncol(Z_gibbs)
  length.alpha <- length(table(clID)) * n.par.rand
  n.par <- length.alpha + ncol(X_gibbs) + 1 + n.par.rand^2
  # MS: Explanation of the number of parameters
  # length.alpha: for each cluster (length(table(clID))) there are n.par.rand random effects variables.
  # ncol(X_gibbs): There are ncol(X_gibbs) fixed effects variables.
  # 1: There is one residual variance.
  # n.par.rand^2: The random effects covariance matrix has n.par.rand * n.par.rand entries
  # (we ignore that the upper triangle is identical to the lower triangle).

  xtx <- solve(t(X_obs) %*% X_obs) #!!! 2016-07-13 war einmal nicht invertierbar!!!
  b.hat.part1 <- xtx %*% t(X_obs)

  data_obs <- data.frame(y_obs = y_obs,
                         clID_obs = clID_obs)
  # MS: get all variables from X and Z into data_obs,
  # but only once (often Z is a subset of X and Zs variables should not appear again in the data set
  # if they are already brought to data_obs by X).

  data_obs[, colnames(X_obs)] <- X_obs
  tmp <- names(Z_obs)[!(names(Z_obs) %in%  colnames(X_obs))]
  data_obs[, tmp] <- Z_obs[, tmp]
  rm(tmp)

  tmp_model_formula <- formula(paste("y_obs ~ 0 + ",
                                     paste(colnames(X_obs), collapse = " + "),
                                     "+ (0 + ",
                                     paste(colnames(Z_obs), collapse = " + "), "| clID_obs)"))
#DEBUG: EINE INTERAKTION BEHEBT DAS PROBLEM AUCH NICHT!!!
 # tmp_model_formula <- formula(y_obs ~ 0 + Intercept + X * cl.var + (0 + Intercept + X | clID_obs))

  #ohne explizite intercept variable änder tauch nichts.
  #tmp_model_formula <- formula(y_obs ~ X * cl.var + (1 + X | clID_obs))
  # MS: Was soll hier passieren???
  res_obs <- lme4::lmer(tmp_model_formula, data = data_obs) #MS: Die Interaktion X * cl.var
  #wird hier nicht beruecksichtig!!!
  vcov.rand <- lme4::VarCorr(res_obs)[[1]][,]
  hyper.df1 <- nclass + dim(vcov.rand)
  hyper.df2 <- vcov.rand * dim(vcov.rand)

  sims <- array(NA, c(n.chains, max.iter, n.par))

  ###generate starting values
#DEBUGCOUNTER <- 0
#DEBUG_COR <- array(dim = max.iter)
  converged <- FALSE
  iter <- 0
  while(!converged & iter < max.iter){
#    DEBUGCOUNTER <- DEBUGCOUNTER + 1
#    set.seed(DEBUGCOUNTER)
    for (k in 1:n.chains){
      if (iter == 0){

        ###generate starting values
        alpha.new <- matrix(rnorm(n = length.alpha), ncol = n.par.rand)
        beta.new <- rnorm(n = ncol(X_obs))
        sigma.y.new <- runif(1)
        sigma.alpha.new <- solve(rWishart(1, df = n.par.rand + 1, Sigma = diag(n.par.rand))[,,1])
      }

      if (iter != 0){

        # MS: use the latest parameters
        alpha.new <- matrix(sims[k, iter, 1:length.alpha], ncol = n.par.rand)
        beta.new <- sims[k, iter, (length.alpha + 1):(length.alpha + ncol(X_obs))]
        sigma.y.new <- sims[k, iter, length.alpha + ncol(X_obs) + 1]
        sigma.alpha.new <- matrix(sims[k, iter,
                (length.alpha + ncol(X_obs) + 2):(length.alpha + ncol(X_obs) + 1 + n.par.rand^2)],
                ncol = n.par.rand)
      }



      draws <- array(NA, dim = c(n.iter, n.par))

      for (s in 1:n.iter){
#        DEBUGCOUNTER <- DEBUGCOUNTER + 1
#        set.seed(DEBUGCOUNTER)
#        save(Z_obs, alpha.new, b.hat.part1, clID_obs, y_obs, sigma.y.new,xtx, file ="updatesigmaalphacoefcrashes.rda")
#        print(DEBUGCOUNTER)#1904 war letzer angezeigter seed
#    DEBUG_COR[DEBUGCOUNTER] <- cor(alpha.new)[1,2]
#    save(DEBUG_COR, file = "DEBUG_COR.rda")
#        print(paste("Correlation of alpha.new:", round(DEBUG_COR[DEBUGCOUNTER], 3)))
        #jetzt wars 1047
#        load("updatesigmaalphacoefcrashes.rda")

        if(FALSE){
          #pdf("CorrelationAlphanew_Over1000runs.pdf")
          jpeg("CorrelationAlphanew_Over1000runs.jpg")
          plot(x = 2:1000, y = DEBUG_COR[2:1000], type = "l")
          dev.off()
        }

        # Updata fixed effects
        beta.new <- update.beta.coef(Z_obs = Z_obs, alpha.new = alpha.new, b.hat.part1 = b.hat.part1,
                                     clID_obs = clID_obs, y_obs = y_obs,
                                     sigma.y.new = sigma.y.new, xtx = xtx)

        # Update the residual variance
        sigma.y.new <- update.sigma.y.coef(y_obs = y_obs, X_obs = X_obs,
                                           beta.new = beta.new, Z_obs = Z_obs,
                                           alpha.new = alpha.new, clID_obs = clID_obs, n.obs = n.obs)

        # Update random effects covariance matrix
        # !!! 2016-07-13  reciprocal condition number
        #save(hyper.df1, hyper.df2, alpha.new, file ="updatesigmaalphacoefcrashes.rda")
        #MACHT PROBLEME!!!
        #Ich denke das Problem ist die (quasi) perfekte Korrelation in alpha.new
        # (oder in hyper.df2)
        sigma.alpha.new <- update.sigma.alpha.coef(hyper.df1 = hyper.df1,
                                                   hyper.df2 = hyper.df2, alpha.new = alpha.new)

        # Update the cluster specific random effects
        #
        alpha.new <- update.alpha.coef(Z_obs = Z_obs, y_obs = y_obs, X_obs = X_obs,
                                       beta.new = beta.new, clID_obs = clID_obs,
                                       sigma.y.new = sigma.y.new,
                                       sigma.alpha.new.inv = solve(sigma.alpha.new))

        # Store the drawn parameters
        draws[s, ] <- c(alpha.new, beta.new, sigma.y.new, sigma.alpha.new)
      }

      sims[k, (iter + 1):(iter + n.iter), ] <- draws
    }

    iter <- iter + n.iter

    ###discard first part of the MCMC chain

    chains <- sims[, ceiling(iter * burn.in):iter, ]

    ###turn into coda usable object

    get.mcmc <- vector("list", n.chains)

    for (k in 1:n.chains){
      get.mcmc[[k]] <- coda::as.mcmc(chains[k,,])
    }

    ###calculate R.hat
    R.hat <- coda::gelman.diag(get.mcmc, autoburnin = FALSE, multivariate = FALSE)
    converged <- max(R.hat$psrf[,1]) < 1.05
  }

  print(paste("converged after", iter, "iterations."))

  return(chains)
}


# Generate documentation with devtools::document()
# Build package with devtools::build() and devtools::build(binary = TRUE) for zips
