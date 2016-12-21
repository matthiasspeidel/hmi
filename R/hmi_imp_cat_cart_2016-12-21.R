#' The function to impute categorical variables
#'
#' The function uses regression trees for imputation. The principle is the following:
#' For each observation it is calculated at wich leave it would end.
#' Then one (randomly selected) observation of the other observations found on this leave,
#' functions as a donor.
#' @param y_imp_multi A Vector with the variable to impute.
#' @param X_imp_multi A data.frame with the fixed effects variables.
#' @return A n x 1 matrix.
imp_cat_cart <- function(y_imp_multi,
                         X_imp_multi){


  categories <- levels(y_imp_multi)

  #Initialising the returning vector
  y_imp <- as.matrix(y_imp_multi, ncol = 1)

  #the missing indactor indicates, which values of y are missing.
  missind <- is.na(y_imp_multi)


  types <- array(dim = ncol(X_imp_multi))
  for(j in 1:length(types)) types[j] <- get_type(X_imp_multi[, j])
  #need_stand <- types == "cont"
  categorical <- types == "categorical"

  #remove categories with more than 10 observations as the model in the current form
  #will cause later numerical probles
  too_many_levels <- colnames(X_imp_multi[, categorical, drop = FALSE])[
    apply(X_imp_multi[, categorical, drop = FALSE], 2, function(x) nlevels(factor(x))) > 10]
  X_imp_multi <- X_imp_multi[, !names(X_imp_multi) %in% too_many_levels, drop = FALSE]


  n <- length(y_imp_multi)
  lmstart <- lm(rnorm(n) ~ 0 +., data = X_imp_multi)

  X_model_matrix_1 <- model.matrix(lmstart)
  xnames_1 <- paste("X", 1:ncol(X_model_matrix_1), sep = "")

  tmp_1 <- data.frame(y = rnorm(n))
  tmp_1[, xnames_1] <- X_model_matrix_1

  reg_1 <- lm(y ~ 0 + . , data = tmp_1)

  blob <- y_imp_multi
  tmp_2 <- data.frame(y = blob)

  xnames_2 <- xnames_1[!is.na(coefficients(reg_1))]
  tmp_2[, xnames_2] <- X_model_matrix_1[, !is.na(coefficients(reg_1)), drop = FALSE]


  # estimate tree
  my_tree <- tree(formula = y ~ 0 + ., data = tmp_2)

  # get predictions on all data
  predict_current<- predict(my_tree, newdata = list(tmp_2)[[1]],
                             type = "where")

  # get predictions on data with missing target variable.
  predict_current_mis <- predict_current[missind]


  #Da der Schleifen-Indikator l1 die (nicht aufsteigende) Nummerierung der Bl?tter
  #durchgeht, brauchen wir counter als durchlaufenden Z?hler.
  counter <- 0

  #Alle Numerierungen der leaves durchgehen, aus denen gezogen werden muss,
  #weil Beobachtungen mit fehlendem y diesem Blatt zugewiesen worden sind.
  for(l1 in as.numeric(names(table(predict_current_mis)))){

    counter <- counter + 1

    # get the absolute frequencies of the categories in the current leave
    # therefore the relative frequencies have to be multiplied with the number of observations in the leave
    # note: in theory you wouldn't need to round. but to avoid numerical problems like
    # 0.1*3*10 == 3, they are rounded to the next integer.
    abs_freq <- round(my_tree$frame$yprob[l1, ] * my_tree$frame$n[l1])


    # generate the donor
    donor <- rep(categories, times = abs_freq)

    # get sampling probabilities
    my_prob <- runif(length(donor))
    # norming it to sum(my_prob) == 1
    my_prob <- my_prob/sum(my_prob)

    # calculate how many observations with a missing Y are found in this leave
    # (needed to know how many new Ys have to be sampled)
    nmiss <- table(predict_current_mis)[counter]

    # select the observations landed in the current leave and have missing values.
    selection <- (predict_current == l1) & missind

    # now sample values from the donor for the selected observations
    y_imp[selection] <- sample(donor, size = nmiss, replace = TRUE,
                               prob = my_prob)

  }

  return(y_imp)
}
