
#' Cycling
#'
#' Function to do one imputation cycle on the given data. The function cycles through
#' every variable sequentially imputing the values, that are NA in the original data set
#' in that current variable. The function determines the type of the variable
#' and calls the suitable imputation function.
#' @param data_before The data.frame with the variables to impute.
#' @param original_data The original data.frame the user passed to \code{hmi}.
#' @param fe A list with the decomposed elements of the \code{model_formula}.
#' @param interaction_names A list with the names of the variables
#' that have been generated as interaction variables
#' @param list_of_types a list where each list element has the name of a variable
#' in the data.frame. The elements have to contain a single character denoting the type of the variable.
#' See \code{get_type} for details about the variable types.
#' With the function \code{list_of_types_maker}, the user can get the framework for this object.
#' In most scenarios this is should not be necessary.
#' One example where it might be necessary is when only two observations
#' of a continuous variable are left - because in this case \code{get_type}
#' interpret is variable to be binary. Wrong is it in no case.
#' @param NA_locator A n x p matrix localizing the missing values in the original
#' dataset. The elements are TRUE if the original data are missing and FALSE if the
#' are observed.
#' @param nitt An integer defining number of MCMC iterations (see \code{MCMCglmm}).
#' @param burnin burnin A numeric value between 0 and 1 for the desired percentage of
#' Gibbs samples that shall be regarded as burnin.
#' @param thin An integer to set the thinning interval range. If thin = 1,
#' every iteration of the Gibbs-sampling chain will be kept. For highly autocorrelated
#' chains, that are only examined by few iterations (say less than 1000),
#' the \code{geweke.diag} might fail to detect convergence. In such cases it is
#' essential to look a chain free from autocorrelation.
#' @param pvalue A numeric between 0 and 1 denoting the threshold of p-values a variable in the imputation
#' model should not exceed. If they do, they are excluded from the imputation model.
#' @param mn An integer defining the minimum number of individuals per cluster.
#' @param k An integer defining the allowed maximum of levels in a factor covariate.
#' @param spike A numeric value saying which value in the semi-continuous data might be the spike.
#' Or a list with with such values and names identical to the variables with spikes
#' (see \code{list_of_spikes_maker} for details).
#' @param rounding_degrees A numeric vector with the presumed rounding degrees. Or a list with rounding degrees,
#' where each list element has the name of a rounded continuous variable. Such a list can be generated
#' using \code{list_of_rounding_degrees_maker(data)}.
#' @param rounding_covariates A list for each rounded continuous variable with a character vector
#' containing the covariate names from the original rounding formula.
#' The transformation takes place in the wrapper function.
#' @return A data.frame where the values, that have a missing value in the original
#' dataset, are imputed.
#' @export
imputationcycle <- function(data_before,
                            original_data,
                            NA_locator,
                            fe, interaction_names,
                            list_of_types,
                            nitt,
                            burnin,
                            thin,
                            pvalue = 0.2,
                            mn,
                            k = Inf,
                            spike = NULL,
                            rounding_degrees = NULL,
                            rounding_covariates){

  if(!is.data.frame(data_before)){
    stop("data_before has to be a data.frame")
  }

  if(!is.matrix(NA_locator)){
    stop("NA_locator has to be a matrix")
  }

  original_ordering <- colnames(data_before)
  missing_rates <- colMeans(NA_locator)

  #Sort data by missing rate
  data_before <- data_before[, names(sort(missing_rates))[names(sort(missing_rates)) %in% names(data_before)], drop = FALSE]

  NA_locator <- NA_locator[, names(sort(missing_rates)), drop = FALSE]

  #update missing rates after the sorting
  missing_rates <- colMeans(NA_locator)

  #get variables with missing values
  incomplete_variables <- names(missing_rates)[missing_rates > 0]

  variables_to_impute <- union(incomplete_variables, names(list_of_types)[list_of_types == "roundedcont" |
                                                                            list_of_types == "interval"])

  #initialize list with the chains of Gibbs-samples
  chains <- list()

  for(l2 in variables_to_impute){#do the imputation cycle

    #restore the original status of the variable which is to be imputed
    data_before[, l2] <- original_data[, l2]

    # get type
    tmp_type <- list_of_types[l2][[1]]

    #update the interaction variables
    if(fe$interactions != ""){
      for(j in 1:length(fe$interactions)){
        #get the product of the interaction variables
        interaction <- apply(data_before[, strsplit(fe$interactions[j], ":")[[1]], drop = FALSE], 1, prod)

        #add the interaction to the dataset
        data_before[, interaction_names[[j]]] <- interaction

        #This step has to be repeated after each imputation of a variable.
      }
    }else{
      interaction_names <- ""
    }

    #generate temporary data set, that contains only the variable to impute and
    #the covariates where no missing values occurs

    # First X and Z are generated
    # X contains the fixed effects variables and interaction variables
    X <- data_before[, union(fe$fixedeffects_varname[fe$fixedeffects_varname %in% names(data_before)],
                             interaction_names[[1]][interaction_names[[1]] %in% names(data_before)]),
                     drop = FALSE]

    # now safe the variables used for modeling random effects
    Z <- data_before[, fe$randomeffects_varname[fe$randomeffects_varname %in% names(data_before)],
                     drop = FALSE]

    tmp_X <- X[, names(X) != l2, drop = FALSE]
    tmp_Z <- Z[, names(Z) != l2, drop = FALSE]

    #now, remove variables from X an Z that currently have NAs
    #(wich should only occur in the first, initial imputation cycle)
    tmp_X <- tmp_X[, stats::complete.cases(t(tmp_X)), drop = FALSE]

    #To preserve the relation of y and x properly,
    #y has to be included as covariate in the imputation model for x.
    #But not only as a fix effect covariate,
    #but also as a random effect covariate in the cases where x
    #is a random effects covariate in the analysis model for y.
    if(l2 %in% fe$randomeffects_varname){
      tmp_Z <- cbind(tmp_Z, X[, fe$target_varname])
    }

    tmp_Z <- tmp_Z[, stats::complete.cases(t(tmp_Z)), drop = FALSE]
    #check number of observed values in each cluster.

    if(fe$clID_varname != ""){

      # Make sure that there are no missing values in the cluster ID:
      data_before[, fe$clID_varname] <- sample_imp(data_before[, fe$clID_varname])

      tab_1 <- table(data_before[!is.na(data_before[, l2]), fe$clID_varname])

      #merge small clusters to a large big cluster
      safety_counter <- 0
      while(any(tab_1 < mn) & safety_counter < 100){
        safety_counter <- safety_counter + 1

        #step 1: determine the smallest and the second smallest cluster.
        smallest_cluster <- tab_1[tab_1 == sort(as.numeric(tab_1))[1]][1]
        tab_1_without_smallest_cluster <- tab_1[names(tab_1) != names(smallest_cluster)]
        second_smallest_cluster <- tab_1_without_smallest_cluster[
          tab_1_without_smallest_cluster ==
            sort(as.numeric(tab_1_without_smallest_cluster))[1]][1]

        #step 2: new cluster name
        new_name <-  names(second_smallest_cluster)

        levels(data_before[, fe$clID_varname])[levels(data_before[, fe$clID_varname]) ==
                                                 names(smallest_cluster)] <- new_name
        levels(data_before[, fe$clID_varname])[levels(data_before[, fe$clID_varname]) ==
                                                 names(second_smallest_cluster)] <- new_name
        #step 3: repeat.

        tab_1 <- table(data_before[!is.na(data_before[, l2]), fe$clID_varname])

      }

    }


    # if too few observations are left for a usefull imputation,
    # the function stops
    if(sum(!is.na(data_before[, l2])) <= 2){
      stop(paste("Too few observations left in variable", l2))
    }

    #If no complete covariates are left
    #(which currently should not happen, as we currently force X to include an intercept variabel)
    #We have to do sample imputation.
    if(ncol(tmp_X) == 0){
      tmp_type <- "sample_imp"
    }

    #check whether a single level or multilevel imputation should be run
    # the standard is basically a multilevel imputation, but if one (or more) of the following
    # conditions is met, a single level imputation is run:
    # 1. there was no cluster ID found by the extract_varnames function
    # 2. the cluster ID cannot be found in the data.
    # 3. the variable to impute is the cluster ID.
    # 4. tmp_Z is currently empty
    # Explanation for 1 and 2: Without an cluster ID in the formula or the data, no multilevel model.
    # Explanation for 3: multilevel models say that the influence of the covariates changes from cluster to cluster
    # so cluster IDs are part of the covariates and not the target variable.
    # Explanation for 4: It can happen that the only random effects variable is going to be imputed,
    # resulting in an empty tmp_Z. But without random effects variables, no random effects model.
    use_single_level <- fe$clID_varname == "" || !(fe$clID_varname %in% names(data_before)) ||
      ncol(tmp_Z) == 0

    if(is.list(spike)){
      spike_tmp <- spike[[l2]]
    }else{
      spike_tmp <- spike
    }

    if(is.list(rounding_degrees)){
      rounding_degrees_tmp <- rounding_degrees[[l2]]
    }else{
      rounding_degrees_tmp <- rounding_degrees
    }

    if(is.null(tmp_type)){
      tmp_type <- list_of_types_maker(data_before[, l2, drop = FALSE],
                                      spike = spike, rounding_degrees = rounding_degrees)[[1]]
    }

    if(tmp_type == "sample_imp"){
      imp <- sample_imp(data_before[, l2])
    }

    if(tmp_type == "binary"){
      if(use_single_level){

        imp <- imp_binary_single(y_imp = data_before[, l2],
                                 X_imp = tmp_X,
                                 pvalue = pvalue,
                                 k = k)

      }else{

        tmp <- imp_binary_multi(y_imp = data_before[, l2],
                                X_imp = tmp_X,
                                Z_imp = tmp_Z,
                                clID = data_before[, fe$clID_varname],
                                nitt = nitt,
                                burnin = burnin,
                                thin = thin,
                                pvalue = pvalue,
                                k = k)

        imp <- tmp$y_ret
        chains[[l2]]$Sol <- tmp$Sol
        chains[[l2]]$VCV <- tmp$VCV
      }
    }

    if(tmp_type == "cont"){

      if(use_single_level){

        imp <- imp_cont_single(y_imp = data_before[, l2],
                               X_imp = tmp_X,
                               pvalue = pvalue,
                               k = k)

      }else{

        tmp <- imp_cont_multi(y_imp = data_before[, l2],
                              X_imp = tmp_X,
                              Z_imp = tmp_Z,
                              clID = data_before[, fe$clID_varname],
                              nitt = nitt,
                              burnin = burnin,
                              thin = thin,
                              pvalue = pvalue,
                              k = k)
        imp <- tmp$y_ret
        chains[[l2]]$Sol <- tmp$Sol
        chains[[l2]]$VCV <- tmp$VCV
      }

    }

    if(tmp_type == "semicont"){
      if(is.list(spike)){
        spike_tmp <- spike[[l2]]
      }else{
        spike_tmp <- spike
      }

      if(use_single_level){

        imp <- imp_semicont_single(y_imp = data_before[, l2],
                                   X_imp = tmp_X,
                                   spike = spike_tmp,
                                   pvalue = pvalue,
                                   k = k)
      }else{

        tmp <- imp_semicont_multi(y_imp = data_before[, l2],
                                  X_imp = tmp_X,
                                  Z_imp = tmp_Z,
                                  clID = data_before[, fe$clID_varname],
                                  spike = spike_tmp,
                                  nitt = nitt,
                                  burnin = burnin,
                                  thin = thin,
                                  pvalue = pvalue,
                                  k = k)

        imp <- tmp$y_ret
        chains[[l2]]$Sol <- tmp$Sol
        chains[[l2]]$VCV <- tmp$VCV

      }

    }


    if(tmp_type == "roundedcont"){
      #yet, we don't have a multilevel imputation for rounded incomes

      #Set up the data.frame including the variable explaining G, the latent rounding tendency
      vars_usable <- rounding_covariates[[l2]]
      if(is.null(vars_usable)){
        vars_usable <- colnames(data_before)
      }
      vars_usable <- vars_usable[vars_usable %in% colnames(data_before)]
      PSI <- data_before[, vars_usable, drop = FALSE]

      #4 different settings are distinguished:
      #Setting 1: rounding_degrees is NULL
      #Setting 2: rounding_degrees is a vector
      #Setting 3a: rounding_degrees is list, with no element for the current variable
      #Setting 3b: rounding_degrees is list, with a specific vector for the current variable

      #Dependent on the setting, the handling of the current variable differs:
      #Setting 1: Rounding degrees are proposed by hmi.
      #Setting 2: The vector is used.
      #Setting 3a: Rounding degrees are proposed by hmi.
      #Setting 3b: The vector is used.


      if(is.null(rounding_degrees)){       ## Setting 1 ##
        rounding_degrees_tmp <- suggest_rounding_degrees(data_before[, l2])

      }else if(is.integer(rounding_degrees)){      ## Setting 2 ##
        rounding_degrees_tmp <- rounding_degrees

      }else if(is.list(rounding_degrees) & is.null(rounding_degrees[[l2]])){      ## Setting 3a ##
        rounding_degrees_tmp <- suggest_rounding_degrees(data_before[, l2])

      }else if(is.list(rounding_degrees) & !is.null(rounding_degrees[[l2]])){      ## Setting 3b ##
        rounding_degrees_tmp <- rounding_degrees[[l2]]
      }


      imp <- imp_roundedcont(y_df = data_before[, l2, drop = FALSE],
                             X_imp = tmp_X,
                             PSI = PSI,
                             pvalue = pvalue,
                             k = k,
                             rounding_degrees = rounding_degrees_tmp)

    }

    if(tmp_type == "interval"){

      #yet, we don't have a multilevel imputation for interval data
      imp <- imp_interval(y_imp = data_before[, l2],
                          X_imp = tmp_X,
                          pvalue = pvalue,
                          k = k)

    }


    if(tmp_type == "count"){

      if(use_single_level){

        tmp <- imp_count_single(y_imp = data_before[, l2],
                                X_imp = tmp_X,
                                nitt = nitt,
                                burnin = burnin,
                                thin = thin,
                                pvalue = pvalue,
                                k = k)

        imp <- tmp$y_ret
        chains[[l2]]$Sol <- tmp$Sol
        chains[[l2]]$VCV <- tmp$VCV


      }else{

        tmp <- imp_count_multi(y_imp = data_before[, l2],
                               X_imp = tmp_X,
                               Z_imp = tmp_Z,
                               clID = data_before[, fe$clID_varname],
                               nitt = nitt,
                               burnin = burnin,
                               thin = thin,
                               pvalue = pvalue,
                               k = k)

        imp <- tmp$y_ret
        chains[[l2]]$Sol <- tmp$Sol
        chains[[l2]]$VCV <- tmp$VCV

      }

    }

    if(tmp_type == "categorical"){
      data_before[, l2] <- original_data[, l2]
      # if the cluster ID is to be imputed, obviously we cannot run a multilevel model.
      if(use_single_level || l2 == fe$clID_varname){

        imp <- imp_cat_single(y_imp = data_before[, l2],
                              X_imp = tmp_X,
                              pvalue = pvalue,
                              k = k)

      }else{

        tmp <- imp_cat_multi(y_imp = as.factor(data_before[, l2]),
                             X_imp = tmp_X,
                             Z_imp = tmp_Z,
                             clID = data_before[, fe$clID_varname],
                             nitt = nitt,
                             burnin = burnin,
                             thin = thin,
                             pvalue = pvalue,
                             k = k)

        imp <- tmp$y_ret
        chains[[l2]]$Sol <- tmp$Sol
        chains[[l2]]$VCV <- tmp$VCV
      }
    }

    if(tmp_type == "ordered_categorical"){

      if(use_single_level){

        imp <- imp_orderedcat_single(y_imp = data_before[, l2],
                                     X_imp = tmp_X,
                                     pvalue = pvalue,
                                     k = k)

      }else{

        tmp <- imp_orderedcat_multi(y_imp = as.factor(data_before[, l2]),
                                    X_imp = tmp_X,
                                    Z_imp = tmp_Z,
                                    clID = data_before[, fe$clID_varname],
                                    nitt = nitt,
                                    burnin = burnin,
                                    thin = thin,
                                    pvalue = pvalue,
                                    k = k)

        imp <- tmp$y_ret
        chains[[l2]]$Sol <- tmp$Sol
        chains[[l2]]$VCV <- tmp$VCV
      }
    }

    if(tmp_type == "intercept"){
      imp <- sample_imp(data_before[, l2])
    }
    data_before[, l2] <- imp
  }
  data_after <- data_before[, original_ordering, drop = FALSE]
  ret <- list(data_after = data_after, chains = chains)
  return(ret)
}
