#' The function for imputation of binary variables.
#'
#' The function is called by the wrapper.
#' @param y_imp_multi A Vector with the variable to impute.
#' @param X_imp_multi A data.frame with the fixed effects variables.
#' @param M An integer defining the number of imputations that should be made.
#' @return A n x M matrix. Each column is one of M imputed y-variables.
imp_count_single <- function(y_imp_multi,
                       X_imp_multi){

  need_stand <- apply(X_imp_multi, 2, get_type) == "cont"
  X_imp_multi_stand <- X_imp_multi
  X_imp_multi_stand[, need_stand] <- scale(X_imp_multi[, need_stand])
  #X_imp_multi %>% mutate_each_(funs(scale), vars = names(need_stand)[need_stand])
  #generate model.matrix (from the class matrix)
  n <- nrow(X_imp_multi_stand)
  X_model_matrix <- model.matrix(rnorm(n) ~ 0 + ., data = X_imp_multi_stand)
  # Remove ` from the variable names
  colnames(X_model_matrix) <- gsub("`", "", colnames(X_model_matrix))

  # -------------- calling the gibbs sampler to get imputation parameters----


  n <- length(y_imp_multi)
  lmstart <- lm(rnorm(n) ~ 0 +., data = X_imp_multi)

  X_model_matrix_1 <- model.matrix(lmstart)
  xnames_1 <- paste("X", 1:ncol(X_model_matrix_1), sep = "")

  tmp_1 <- data.frame(y = rnorm(n))
  tmp_1[, xnames_1] <- X_model_matrix_1

  reg_1 <- lm(y ~ 0 + . , data = tmp_1)

  blob <- y_imp_multi
  tmp_2 <- data.frame(target = blob)

  xnames_2 <- xnames_1[!is.na(coefficients(reg_1))]
  X_model_matrix_2 <- X_model_matrix_1[, !is.na(coefficients(reg_1)), drop = FALSE]
  tmp_2[, xnames_2] <- X_model_matrix_2

  tmp_2[, "ClID"] <- clID

  fixformula <- formula(paste("target~", paste(xnames_2, collapse = "+"), "- 1", sep = ""))

  prior <- list(R = list(V = 1e-07, nu = -2))

  MCMCglmm_draws <- MCMCglmm::MCMCglmm(fixformula, data = tmp_2,
                                       verbose = FALSE, pr = TRUE, prior = prior,
                                       family = "poisson",
                                       saveX = TRUE,
                                       nitt = 3000,
                                       thin = 10,
                                       burnin = 1000)

  pointdraws <- MCMCglmm_draws$Sol
  xdraws <- pointdraws[, 1:ncol(X_model_matrix_2), drop = FALSE]
  variancedraws <- MCMCglmm_draws$VCV
  # the last column contains the variance (not standard deviation) of the residuals

  number_of_draws <- nrow(pointdraws)
  select.record <- sample(1:number_of_draws, M, replace = TRUE)

  # -------------------- drawing samples with the parameters from the gibbs sampler --------
  y_imp <- array(NA, dim = c(n, M))
  ###start imputation
  for (j in 1:M){

    fix.eff.imp <- matrix(xdraws[select.record[j], ], nrow = ncol(X_model_matrix_2))

    sigma.y.imp <- sqrt(variancedraws[select.record[j], ncol(variancedraws)])


    lambda <- exp(rnorm(n, X_model_matrix_2 %*% fix.eff.imp, sigma.y.imp))

    #problematic <- lambda < 0
    #while(any(problematic)){
    #  lambda[problematic] <- rnorm(sum(problematic),
    #                               X_model_matrix_2[problematic, , drop = FALSE] %*% fix.eff.imp +
    #                                 apply(Z_imp_multi_stand[problematic, , drop = FALSE] * rand.eff.imp[clID,], 1, sum),
    #                               sigma.y.imp)
    #  problematic <- lambda < 0
    #}


    y_imp[, j] <- ifelse(is.na(y_imp_multi), rpois(n, lambda), y_imp_multi)
  }

  # --------- returning the imputed data --------------
  return(y_imp)

  #end testfield
}



# Generate documentation with devtools::document()
# Build package with devtools::build() and devtools::build(binary = TRUE) for zips
