# get fixed effect predictions from an lmer model

get_prediction_frame <- function(data, linear_predictors = NULL
                                 , fixed_predictors = NULL
                            , factor_predictors = NULL){
  
  # linear_predictors and fixed_predictors are character vectors of the 
  # variables which want to be represented as linear(11 values from max to min)
  # , and the mean values, respectively.
  # factor_predictors should be a named list with the unique values of the
  # factor to plot
  
  if(!require(rlist)){library(rlist)}
  
  # if(is.null(linear_predictors)){
  #   linear_predictors <- c(1:2)
  # } # Make a dummy variable to make life easy. Code will run but nothing 
  # # will happen with this output in the end
  
  # if(is.null(fixed_predictors)){
  #   fixed_predictors <- c(1:2)
  # } # Make a dummy variable to make life easy. Code will run but nothing 
  # # will happen with this output in the end
  
  # Get the mean value of all fixed (numeric) variables
  fixed_pred_values <- lapply(fixed_predictors, FUN = function(x){
    with(data, mean(get(x)))
  })
  
  # Get a sequence of 100 values from min to max of linear predictors
  linear_pred_values <- lapply(linear_predictors
                               , FUN = function(x){
                                 with(data
                                      , seq(min(get(x))
                                            , max(get(x))
                                            , length.out = 100)
                                 )
                               })
  
  # Make a named list of these values, and pass to expand.grid
  predvals <- list.flatten(list(linear_pred_values, fixed_pred_values))
  if(!is.null(factor_predictors)){predvals <- list.flatten(list(predvals, factor_predictors))}
   if(!is.null(factor_predictors)){
    columnames <- c(linear_predictors, fixed_predictors
      , names(factor_predictors))
  } else {
    columnames <- c(linear_predictors, fixed_predictors)
  }
  
  names(predvals) <- c(linear_predictors, fixed_predictors
                       , names(factor_predictors))
  preddata <- expand.grid(predvals)
  
  # Return the output
  return(preddata)
  
}

# Handy function for getting confidence and prediction intervals
get_fixed_predictor_formula <- function(mod){
  as.formula(paste0("~ ", sub(".*~", "", lme4::nobars(formula(mod)))[3]))
}
