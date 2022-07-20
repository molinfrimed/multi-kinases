library(parallel)
library(naivebayes)
library(caret)
library(e1071)
library(C50)
library(randomForest)
library(neuralnet)
library(xgboost)
library(class)
library(rpart)
library(rpart.plot)
library(glmnet)
library(ROCR)

get.auc <- function(pred){
  auc <- performance(pred, measure = "auc")
  return(unlist(auc@y.values))
}

save.roc <- function(pred, img.path=NULL, size=4, add=FALSE, col="blue"){
  roc <- performance(pred,"tpr","fpr")
  
  if(!is.null(img.path) && !add){
    png(filename = img.path, width = size, height = size, units= "in", res=300)
    plot(roc, add=add, col=col, lwd = 1.5)
    abline(a=0,b=1, lty=2)
    dev.off()
  } else{
    plot(roc, add=add, col=col, lwd = 1.5)
    abline(a=0,b=1, lty=2)
  }
  
}

# rangeNormalization <- function(ds){
#   norm.ds <- apply(ds, 2, function(feat){
#     m <- min(feat, na.rm = TRUE)
#     M <- max(feat, na.rm = TRUE)
#     
#     norm.feat <- (feat - m) / (M - m)
#     return(norm.feat)
#   })
#   return(norm.ds)
# }

get.partition <- function(label, p = 0.7){
  values <- unique(descrs.all$LABEL)
  
  res <- lapply(values, function(val){
    oids <- which(label == val)
    v <- label[oids]
    map.table <- data.frame(orig=oids, new=1:length(v))
    part.ids <- createDataPartition(v, p = p, list = FALSE)
    return(map.table[part.ids, 1])
  })
  
  return(unlist(res))
}

get.descr.by.prot <- function(descr.ds, prot){
  return(descr.ds[which(descr.ds$Uniprot_ID == prot), ])
}

get.boruta.feat <- function(comp.ds, prot.id){
  ids <- which(comp.ds$Uniprot_id == prot.id)
  bonly <- which(comp.ds[ids, "Boruta"])
  descrs <- comp.ds[bonly, "Descriptor"]
  return(descrs)
}

get.lm.feat <- function(comp.ds, prot.id){
  ids <- which(comp.ds$Uniprot_id == prot.id)
  lmonly <- which(comp.ds[ids, "Lm"])
  descrs <- comp.ds[lmonly, "Descriptor"]
  return(descrs)
}

get.common.feat <- function(comp.ds, prot.id){
  ids <- which(comp.ds$Uniprot_id == prot.id)
  commonly <- which(comp.ds[ids, "Label"] == "COMMON")
  descrs <- comp.ds[commonly, "Descriptor"]
  return(descrs)
}

get.union.feat <- function(comp.ds, prot.id){
  ids <- which(comp.ds$Uniprot_id == prot.id)
  descrs <- comp.ds[ids, "Descriptor"]
  return(descrs)
}

# naive bayes naivebayes
naiveModel <- function(train.set, weights=FALSE){
  tryCatch({
    # compute weights for imbalanced dataset
    if(weights){
      nr <- nrow(train.set)
      sumwpos <- sum(train.set$LABEL == "ACTIVE")/nr
      sumwneg <- sum(train.set$LABEL == "INACTIVE")/nr
      model_weights <- c(sumwpos, sumwneg)
    } else{
      model_weights <- NULL
    }
    
    # launch model
    return(naive_bayes(LABEL ~ ., train.set, usekernel = TRUE, prior = model_weights))
  }, error = function(e) {
    message(paste("NAIVE model raised an error for protein", prot))
    message(e)
    return(NULL)
  })
}

# logistic regression e1071
logisticModel <- function(train.set, weights=FALSE){
  tryCatch({
    if(weights){
      model_weights <- ifelse(train.set$LABEL == "ACTIVE",
                              (1/sum(train.set$LABEL == "ACTIVE")) * 0.5,
                              (1/sum(train.set$LABEL == "INACTIVE")) * 0.5)
    } else{
      model_weights <- NULL
    }
    
    # launch model
    return(glm(LABEL ~ ., data = train.set, family = binomial, weight = model_weights))
  }, error = function(e) {
    message(paste("LOGISTIC model raised an error for protein", prot))
    message(e)
    return(NULL)
  })
}

# SVM e1071
svmModel <- function(train.set, test.set, probability=FALSE, weights=FALSE){
  # tuning.params.linear <- tune(svm, LABEL ~., data = train.set,
  #                            ranges = list(cost = seq(1e-4, 10, length.out=5),
  #                                          epsilon = seq(0.001, 0.1, length.out=5)), # gamma coef0 degree epsilon
  #                            tunecontrol = tune.control(sampling = "fix", performances=TRUE),
  #                            validation.x = test.set[, -which(colnames(test.set) == "LABEL")],
  #                            validation.y = test.set$LABEL,
  #                            kernel = "linear"
  # )
  # 
  # svm.model.linear <- svm(LABEL ~., train.set, scale = FALSE,
  #                       cost = tuning.params.linear$best.parameters$cost,
  #                       epsilon = tuning.params.linear$best.parameters$epsilon,
  #                       kernel = "linear")
  # 
  # tuning.params.poly <- tune(svm, LABEL ~., data = train.set,
  #       ranges = list(gamma = seq(1e-3,2, length.out=5),
  #                     coef0 = seq(0, 10, length.out=5),
  #                     degree = c(2:5),
  #                     cost = seq(1e-4, 10, length.out=5),
  #                     epsilon = seq(0.001, 0.1, length.out=5)), # gamma coef0 degree epsilon
  #       tunecontrol = tune.control(sampling = "fix", performances=TRUE),
  #       validation.x = test.set[, -which(colnames(test.set) == "LABEL")],
  #       validation.y = test.set$LABEL,
  #       kernel = "polynomial"
  # )
  # 
  # svm.model.poly <- svm(LABEL ~., train.set, scale = FALSE,
  #                  gamma = tuning.params.poly$best.parameters$gamma,
  #                  cost = tuning.params.poly$best.parameters$cost,
  #                  coef0 = tuning.params.poly$best.parameters$coef0,
  #                  degree = tuning.params.poly$best.parameters$degree,
  #                  epsilon = tuning.params.poly$best.parameters$epsilon,
  #                  kernel = "polynomial")
  # tuning.params.sigmo <- tune(svm, LABEL ~., data = train.set,
  #      ranges = list(gamma = seq(1e-3,2, length.out=5),
  #                    coef0 = seq(0, 10, length.out=5),
  #                    cost = seq(1e-4, 10, length.out=5),
  #                    epsilon = seq(0.001, 0.1, length.out=5)),
  #      tunecontrol = tune.control(sampling = "fix", performances=TRUE),
  #      validation.x = test.set[, -which(colnames(test.set) == "LABEL")],
  #      validation.y = test.set$LABEL,
  #      kernel = "sigmoid"
  #                                           
  # ) # gamma coef0 cost epsilon
  # 
  # svm.model.sigmo <- svm(LABEL ~., train.set, scale = FALSE,
  #                  gamma = tuning.params.sigmo$best.parameters$gamma,
  #                  coef0 = tuning.params.sigmo$best.parameters$coef0,
  #                  cost = tuning.params.sigmo$best.parameters$cost,
  #                  epsilon = tuning.params.sigmo$best.parameters$epsilon,
  #                  kernel = "sigmoid")
  tryCatch({
    # compute weights for imbalanced dataset
    if(weights){
      nr <- nrow(train.set)
      sumwpos <- sum(train.set$LABEL == "ACTIVE")/nr
      sumwneg <- sum(train.set$LABEL == "INACTIVE")/nr
      model_weights <- c(ACTIVE=sumwpos, INACTIVE=sumwneg)
    } else{
      model_weights <- NULL
    }
    
    # launch model
    tuning.params.radial <- tune(svm, LABEL ~., data = train.set,
                                 ranges = list(gamma = seq(1e-3,2, length.out=5), 
                                               cost = seq(1e-4, 10, length.out=5),
                                               epsilon = seq(0.001, 0.1, length.out=5)),
                                 tunecontrol = tune.control(sampling = "fix", performances=TRUE),
                                 validation.x = test.set[, -which(colnames(test.set) == "LABEL")],
                                 validation.y = test.set$LABEL,
                                 kernel = "radial", scale = FALSE) # gamma, cost, epsilon
    svm.model.radial <- svm(LABEL ~., train.set, scale = FALSE,
                            gamma = tuning.params.radial$best.parameters$gamma,
                            cost = tuning.params.radial$best.parameters$cost,
                            epsilon = tuning.params.radial$best.parameters$epsilon,
                            kernel = "radial",
                            probability = probability,
                            class.weights = model_weights)
    return(svm.model.radial)
  }, error = function(e) {
    message(paste("SVM model raised an error for protein", prot))
    message(e)
    return(NULL)
  })
}

# svmModel.prob <- function(train.set, test.set, weights=FALSE){
#   tryCatch({
#     tuning.params.radial <- tune(svm, LABEL ~., data = train.set,
#                                  ranges = list(gamma = seq(1e-3,2, length.out=5), 
#                                                cost = seq(1e-4, 10, length.out=5),
#                                                epsilon = seq(0.001, 0.1, length.out=5)),
#                                  tunecontrol = tune.control(sampling = "fix", performances=TRUE),
#                                  validation.x = test.set[, -which(colnames(test.set) == "LABEL")],
#                                  validation.y = test.set$LABEL,
#                                  kernel = "radial", scale = FALSE) # gamma, cost, epsilon
#     svm.model.radial <- svm(LABEL ~., train.set, scale = FALSE,
#                             gamma = tuning.params.radial$best.parameters$gamma,
#                             cost = tuning.params.radial$best.parameters$cost,
#                             epsilon = tuning.params.radial$best.parameters$epsilon,
#                             kernel = "radial",
#                             probability=TRUE)
#     return(svm.model.radial)
#   }, error = function(e) {
#     message(paste("SVM model raised an error for protein", prot))
#     message(e)
#     return(NULL)
#   })
# }

# decision tree C50
c50Model <- function(train.set, test.set, weights=FALSE){
  tryCatch({
    # compute weights for imbalanced dataset
    if(weights){
      model_weights <- ifelse(train.set$LABEL == "ACTIVE",
                              (1/sum(train.set$LABEL == "ACTIVE")) * 0.5,
                              (1/sum(train.set$LABEL == "INACTIVE")) * 0.5)
    } else {
      model_weights <- NULL
    }
    tuning.params <- tune(C5.0, LABEL ~., data = train.set, 
                          ranges = list(trials = seq(1, 50, length.out=5)),
                          tunecontrol = tune.control(sampling = "fix", performances=TRUE),
                          validation.x = test.set[, -which(colnames(test.set) == "LABEL")],
                          validation.y = test.set$LABEL)
    c50.model <- C5.0(x = train.set[, -which(colnames(train.set) == "LABEL")], y = train.set$LABEL,
                      trials = tuning.params$best.parameters$trials, rules = TRUE,
                      weights = model_weights)
    return(c50.model)
  }, error = function(e) {
    message(paste("DECISION TREE model raised an error for protein", prot))
    message(e)
    return(NULL)
  })
}

# random forest randomForest tunRF
rfModel <- function(train.set, weights=FALSE){
  tryCatch({
    # compute weights for imbalanced dataset
    if(weights){
      nr <- nrow(train.set)
      sumwpos <- sum(train.set$LABEL == "ACTIVE")/nr
      sumwneg <- sum(train.set$LABEL == "INACTIVE")/nr
      model_weights <- c(sumwpos, sumwneg)
    } else{
      model_weights <- NULL
    }
    
    # launch model
    rf.model <- randomForest(LABEL ~ ., data = train.set, classwt = model_weights)
    return(rf.model)
  }, error = function(e) {
    message(paste("RANDOM FOREST model raised an error for protein", prot))
    message(e)
    return(NULL)
  })
}

# neural network neuralnet
neuralModel <- function(train.set, weights=FALSE){
  tryCatch({
    
    # compute weights for imbalanced dataset
    if(weights){
      model_weights <- ifelse(train.set$LABEL == "ACTIVE",
                              (1/sum(train.set$LABEL == "ACTIVE")) * 0.5,
                              (1/sum(train.set$LABEL == "INACTIVE")) * 0.5)
    } else{
      model_weights <- NULL
    }
    
    # launch model and perform tune grid
    set.seed(128)
    # ctrl <- trainControl(method="repeatedcv", repeats = 10)
    # nnetGrid <-  expand.grid(size = 1:5, 
    #                         decay = seq(0.1, 0.5, by = 0.1))
    # neural.model <- train(LABEL ~ ., data = train.set, method = "nnet",
    #       trControl = ctrl, weights = model_weights,
    #       tuneGrid = nnetGrid, metric = "Accuracy")
    neural.model <- neuralnet(LABEL ~ ., data = train.set, hidden = c(3,2))
    return(neural.model)
  }, error = function(e) {
    message(paste("NEURAL NETWORK model raised an error for protein", prot))
    message(e)
    return(NULL)
  })
}

# xgboost
xgboostModel <- function(train.set, weights=FALSE){
  tryCatch({
    # compute weights for imbalanced dataset
    if(weights){
      nr <- nrow(train.set)
      sumwpos <- sum(train.set$LABEL == "ACTIVE")/nr
      sumwneg <- sum(train.set$LABEL == "INACTIVE")/nr
      params <- list("scale_pos_weight" = sumwneg / sumwpos)
    } else{
      params <- list()
    }
    
    # launch model
    id <- which(colnames(train.set) == "LABEL")
    label.set <- ifelse(train.set$LABEL == "ACTIVE", 1, 0)
    xgboost.cv <- xgb.cv(data = data.matrix(train.set[,-id]), label = label.set,
                         nrounds = 500, max_depth = 5, eta = 0.1, nfold = 5)
    totrounds <- which.min(xgboost.cv$evaluation_log$train_rmse_mean)
    xgboost.model <- xgboost(params = params, data = data.matrix(train.set[,-id]), label = label.set,
                             nrounds = totrounds, max_depth = 5, eta = 0.1,
                             objective = "binary:logistic")
    return(xgboost.model)
  }, error = function(e) {
    message(paste("XGBOOST model raised an error for protein", prot))
    message(e)
    return(NULL)
  })
}

# knn class
knnModel <- function(train.set, test.set, k=5, weights=FALSE){
  tryCatch({
    # compute weights for imbalanced dataset
    if(weights){
      model_weights <- ifelse(train.set$LABEL == "ACTIVE",
                              (1/sum(train.set$LABEL == "ACTIVE")) * 0.5,
                              (1/sum(train.set$LABEL == "INACTIVE")) * 0.5)
    } else{
      model_weights <- NULL
    }
    
    # launch model
    ctrl <- trainControl(method="repeatedcv", repeats = 10)
    # knnGrid <-  expand.grid(kmax = 1:50,
    #                         distance = 1:20 ,
    #                         kernel = c('gaussian',  # different weighting types in kknn
    #                                    'triangular',
    #                                    'rectangular',
    #                                    'epanechnikov',
    #                                    'optimal'))
    knn.model <- train(LABEL ~ ., data = train.set, method = "knn",
                       trControl = ctrl, weights = model_weights
                       #,tuneGrid = knnGrid
                       )
    return(knn.model)
  }, error = function(e) {
    message(paste("KNN model raised an error for protein", prot))
    message(e)
    return(NULL)
  })
}

# cart rpart
cartModel <- function(train.set, test.set, weights=FALSE){
  tryCatch({
    if(weights){
      nr <- nrow(train.set)
      sumwpos <- sum(train.set$LABEL == "ACTIVE")/nr
      sumwneg <- sum(train.set$LABEL == "INACTIVE")/nr
      parms = list(prior=c(sumwpos,sumwneg))
    } else{
      parms <- NULL
    }
    
    # launch model
    minlabel <- min(table(test.set$LABEL))
    best.cart.model <- tune.rpart(LABEL ~ ., data = train.set,
                                  minsplit = seq(5, minlabel, length.out=5),
                                  cp = seq(0.001, 0.01, length.out=5),
                                  maxcompete = c(1,2,3,4,5),
                                  maxdepth = seq(1,30, length.out=5))
    cart.model <- rpart(LABEL ~ ., data = train.set,
                        minsplit = best.cart.model$best.parameters$minsplit,
                        cp = best.cart.model$best.parameters$cp,
                        maxdepth = best.cart.model$best.parameters$maxdepth,
                        maxcompete = best.cart.model$best.parameters$maxcompete,
                        parms = parms)
    return(cart.model)
  }, error = function(e) {
    message(paste("CART model raised an error for protein", prot))
    message(e)
    return(NULL)
  })
}

# lasso regression glmnet
lassoModel <- function(train.set, test.set, weights=FALSE){
  tryCatch({
    if(weights){
      model_weights <- ifelse(train.set$LABEL == "ACTIVE",
                              (1/sum(train.set$LABEL == "ACTIVE")) * 0.5,
                              (1/sum(train.set$LABEL == "INACTIVE")) * 0.5)
    } else{
      model_weights <- NULL
    }
    
    # launch model
    id <- which(colnames(train.set) == "LABEL")
    label.set <- ifelse(train.set$LABEL == "ACTIVE", 1, 0)
    cv.res <- cv.glmnet(as.matrix(train.set[,-id]), label.set)
    label.test <- ifelse(test.set$LABEL == "ACTIVE", 1, 0)
    best.lasso.model <- glmnet(as.matrix(train.set[,-id]), label.set, alpha = 1,
                               family = "binomial", standardize = FALSE,
                               lambda = cv.res$lambda.min,
                               weights = model_weights)
    return(best.lasso.model)
  }, error = function(e) {
    message(paste("LASSO model raised an error for protein", prot))
    message(e)
    return(NULL)
  })
}

# ridge regression glmnet
ridgeModel <- function(train.set, test.set, weights=FALSE){
  tryCatch({
    if(weights){
      model_weights <- ifelse(train.set$LABEL == "ACTIVE",
                              (1/sum(train.set$LABEL == "ACTIVE")) * 0.5,
                              (1/sum(train.set$LABEL == "INACTIVE")) * 0.5)
    } else{
      model_weights <- NULL
    }
    
    # launch model
    id <- which(colnames(train.set) == "LABEL")
    label.set <- ifelse(train.set$LABEL == "ACTIVE", 1, 0)
    cv.res <- cv.glmnet(as.matrix(train.set[,-id]), label.set)
    label.test <- ifelse(test.set$LABEL == "ACTIVE", 1, 0)
    best.ridge.model <- glmnet(as.matrix(train.set[,-id]), label.set, alpha = 0,
                               family = "binomial", standardize = FALSE,
                               lambda = cv.res$lambda.min,
                               weights = model_weights)
    return(best.ridge.model)
  }, error = function(e) {
    message(paste("RIDGE model raised an error for protein", prot))
    message(e)
    return(NULL)
  })
}

# elasticnet regression glmnet
elnetModel <- function(train.set, weights=FALSE){
  tryCatch({
    # compute weights for imbalanced dataset
    if(weights){
      model_weights <- ifelse(train.set$LABEL == "ACTIVE",
                              (1/sum(train.set$LABEL == "ACTIVE")) * 0.5,
                              (1/sum(train.set$LABEL == "INACTIVE")) * 0.5)
    } else {
      model_weights <- NULL
    }
    
    # launch model
    elnet.model <- train(LABEL ~., data = train.set, method = "glmnet",
                         family="binomial",
                         trControl = trainControl("cv", number = 10),
                         tuneLength = 10, weights = model_weights)
    return(elnet.model)
  }, error = function(e) {
    message(paste("ELASTICNET model raised an error for protein", prot))
    message(e)
    return(NULL)
  })
}

nullPerformances <- function(){
  return(list(overall=c(
    Accuracy=NA, Kappa=NA, AccuracyLower=NA, AccuracyUpper=NA, AccuracyNull=NA,
    AccuracyPValue=NA, McnemarPValue=NA),
    byClass = c(
      Sensitivity=NA, Specificity=NA, "Pos Pred Value"=NA, "Neg Pred Value"=NA, 
      Precision=NA, Recall=NA, F1=NA, Prevalence=NA,
      "Detection Rate"=NA, "Detection Prevalence"=NA, "Balanced Accuracy"=NA))
  )
}
