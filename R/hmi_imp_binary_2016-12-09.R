#' The function for imputation of binary variables.
#'
#' The function is called by the wrapper.
#' @param y_imp_multi A Vector with the variable to impute.
#' @param X_imp_multi A data.frame with the fixed effects variables.
#' @param M An integer defining the number of imputations that should be made.
#' @return A n x M matrix. Each column is one of M imputed y-variables.
imp_binary <- function(y_imp_multi,
                      X_imp_multi,
                      M = 10){

  # -----------------------------preparing the data ------------------
  # -- standardise the covariates in X (which are numeric and no intercept)

  types <- array(dim = ncol(X_imp_multi))
  for(j in 1:length(types)) types[j] <- get_type(X_imp_multi[, j])

  need_stand <- types == "cont" # apply(X_imp_multi, 2, get_type) didnt work
  #STANDARDISING SEMICONT AND ROUNDEDCONT, WOULD NEED TO CHANGE THEIR IMPUTATION ROUTINES!
  X_imp_multi_stand <- X_imp_multi
  #X_imp_multi_stand[, need_stand] <- scale(X_imp_multi[, need_stand])
  #X_imp_multi %>% mutate_each_(funs(scale), vars = names(need_stand)[need_stand])
  #generate model.matrix (from the class matrix)
  n <- nrow(X_imp_multi_stand)
  tmp0 <- data.frame(y_binary = rnorm(n), X_imp_multi_stand)

  X_model_matrix_1 <- model.matrix(y_binary ~ 0 +., data = tmp0)
  #PROBLEM If there is an intercept variable present!!!
  #explanation: categorical with k Kategories make (normally) only k-1 dummy variables in the data
  #but not if the formula says (0+...), then k dummy variables are build.
  #if then an additional intercept variable is present, the model will produce an NA estimator for it
  #SO I NEED AN ADDITION STEP WHICH EXCLUDES NA VARIABLES!!! OR I PREVENT NAs IN FOREHAND!!!
  #!!! BESSER LOESEN!!!
    #alt: model.matrix(rnorm(n) ~ 0 + ., data = X_imp_multi_stand)
  # Remove ` from the variable names
  colnames(X_model_matrix_1) <- gsub("`", "", colnames(X_model_matrix_1))


  # -------------- calling the gibbs sampler to get imputation parameters----

  tmp_1 <- data.frame(y_binary = y_imp_multi)
  xnames_1 <- paste("X", 1:ncol(X_model_matrix_1), sep = "")
  tmp_1[, xnames_1] <- X_model_matrix_1

  fixformula_1 <- formula(paste("y_binary ~ 0 +", paste(xnames_1, collapse = "+"), sep = ""))

  reg_1 <- glm(fixformula_1, data = tmp_1, family = binomial(link = "logit"))

  # remove those variables where no effect can be estimated.
  # e.g. if an intercept variable is present, but the model is y~ 0 + Intercept,
  # lm / glm include a dummy variable for every category of a categorical variable (and not k-1).
  # in other words: lm/glm doesn't realise "Intercept" as the intercept variable.
  tmp_2 <- data.frame(y_binary = y_imp_multi)

  xnames_2 <- xnames_1[!is.na(coefficients(reg_1))]
  tmp_2[, xnames_2] <- X_model_matrix_1[, !is.na(coefficients(reg_1)), drop = FALSE]

  fixformula_2 <- formula(paste("y_binary ~ 0 +", paste(xnames_2, collapse = "+"), sep = ""))
  reg_2 <- glm(fixformula_2, data = tmp_2, family = binomial(link = "logit"))
  X_model_matrix_2 <- model.matrix(reg_2)

  reg_coef <- coefficients(reg_2)
  reg_coef_var <- vcov(reg_2)


  linkfunction <- function(x){
    ret <- boot::inv.logit(x)
    return(ret)
  }

  y_imp <- array(NA, dim = c(n, M))
  ###start imputation
  for (j in 1:M){

    fix_eff_imp <- matrix(newbeta <- MASS::mvrnorm(1, reg_coef, reg_coef_var),
                          nrow = ncol(X_model_matrix_2))

    linearpredictor <- rnorm(n, X_model_matrix_2 %*% fix_eff_imp, 0)


    one_prob <- linkfunction(linearpredictor)


    y_temp <- as.numeric(runif(n) < one_prob)

    y_imp[, j] <- ifelse(is.na(y_imp_multi), y_temp, y_imp_multi)
  }

  # --------- returning the imputed data --------------
  return(y_imp)
}


# Generate documentation with devtools::document()
# Build package with devtools::build() and devtools::build(binary = TRUE) for zips
