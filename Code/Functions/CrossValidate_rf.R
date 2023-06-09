#### Make a function for Random Forest cross validation              ####
####     --find optimal parameters and get cross-validated R-squared ####

CrossValidate_rf <- function(y, data
                             , k = NULL
                             , repeats = NULL
                             , ntree = NULL
                             , mtry = NULL){
  
  ## (1) Set things up
  
  # Install and load caret and randomForest packages
  if (!require(caret)) install.packages('caret')
  library(caret)
  if (!require(randomForest)) install.packages('randomForest')
  library(randomForest)
  
  if(is.null(k)){k <- 10} # Number of data subsets to try with
  if(is.null(repeats)){repeats <- 10} # Number of repeats for each subset 
  if(is.null(ntree)){ntree <- 1000} # Number of decision trees in the random forest
  if(is.null(mtry)){mtry <- 2:(ncol(data)-1)}
  
  # Create a formula
  RF_formula <- as.formula(paste(y, "~ ."))
  
  # Number of predictors to try
  tune_grid <- expand.grid(.mtry=mtry)
  
  # Set up training sets - 3 fold cross validation (n for each sample = 4) 
  train_control <- trainControl(method = "repeatedcv"
                                , number = k # How many folds -- given n = 14 this has to be small
                                , repeats = repeats # # Do each fold how many times?
  )
  
  ## (2) Fit the random forest model
  model <- train(RF_formula
                 , data = data
                 , trControl = train_control
                 , tuneGrid = tune_grid
                 , method = "rf"
  )
  
  ## (3) Return results of the cross validation and model optimization
  
  # Easy to read cross validation results
  print(model)
  FullResults <- model$results[, c(1:3, 6)]
  
  # Plot the effect of number of predictors on model performance
  plot(Rsquared ~ mtry, data = FullResults
       , pch = 16
       , ylim = c(0,1))
  
  # Error Bars
  with(FullResults,
       arrows(mtry, Rsquared-RsquaredSD
              , mtry, Rsquared+RsquaredSD
              , length = 0.05
              , angle = 90
              , code = 3)
  )
  
  ModelPerformance <- recordPlot()
  
  ## (4) Take the best model as found by cross validation
  rf_best <- randomForest(RF_formula
                          , data = data
                          , ntree = 1000
                          , mtry = model$bestTune$mtry
  )
  
  ## (5) Plot and print the results
  varImpPlot(rf_best, main = y)
  ImportancePlot <- recordPlot()
  
  ImportanceMetric <- importance(rf_best)
  ImportanceMetric <- data.frame(Order = row.names(ImportanceMetric)[order(ImportanceMetric)]
                                 , IncNodePurity = ImportanceMetric[order(ImportanceMetric)])
  
  return(list(rf_best = rf_best
              , CrossValidation_Results = FullResults
              , mtry_Plot = ModelPerformance
              , Variable_Importance = ImportanceMetric))
  
}

# Byte compile this to make it faster
CrossValidate_rf <- compiler::cmpfun(CrossValidate_rf)
