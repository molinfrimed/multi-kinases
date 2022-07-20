options(stringsAsFactors = FALSE)

### Define path variables
SCRIPT.DIR <- "scripts"
DATA.DIR <- "data"
MODELS.DIR <- "../Version_2/RESULTS/"

# Load functions and scaled data (molecular descriptors computed with PaDEL)
# Example code:
#   descriptors <- data.frame(scale(raw.descriptors[,-1]))
#   rownames(descriptors) <- raw.descriptors$Name
# 
#   # replace NaN values with 0
#   descriptors[,names(which(apply(descriptors, 2, function(col){
#     all(is.nan(col))
#   })))] <- 0
source(file.path(SCRIPT.DIR, "model_functions.R"))
load(file.path(DATA.DIR, "ligands_descriptors.RData"))

# Load best models and set metric
best.models <- read.delim(file = file.path(DATA.DIR, "best_mfs_kuala.tsv"))
metric <- "F1"
res.dir <- "kuala_results"

# Create complete set of ligands
proteins <- unique(x = best.models$prot)

# Load feature set comparison
final.comp <- read.delim(file = file.path(DATA.DIR,"feat_sel_df.tsv"))

# Create results folder in the current working directory
dir.create(file.path(res.dir), showWarnings = FALSE)

# Predict activity and collect the results
res <- lapply(proteins, function(prot){
  cat("Selected ", prot, "\n")
  
  # Select best model and feature set
  mid <- which(best.models$prot == prot & best.models$metric == metric)
  featset.name <- best.models$featset[mid]
  bmodel <- tolower(best.models$method[mid])
  
  cat("Best model for", prot, ":", bmodel, "\n",
      file = paste0(res.dir, "/log.txt"), append = TRUE)
  
  # Set feature set from user selection
  if(featset.name == "boruta"){
    # select boruta descriptor
    boruta.feat <- get.boruta.feat(final.comp, prot)
    feature.set <- descriptors[, boruta.feat]
    cat("  Feature set name:", featset.name, ", number:", ncol(feature.set), "\n", 
        file = paste0(res.dir, "/log.txt"), append = TRUE)
  } else if(featset.name == "linear"){
    # select lm descriptor
    lm.feat <- get.lm.feat(final.comp, prot)
    feature.set <- descriptors[, lm.feat]
    cat("  Feature set name:", featset.name, ", number:", ncol(feature.set), "\n",
        file = paste0(res.dir, "/log.txt"), append = TRUE)
  } else if(featset.name == "common"){
    # select common descriptor
    common.feat <- get.common.feat(final.comp, prot)
    feature.set <- descriptors[, common.feat]
    cat("  Feature set name:", featset.name, ", number:", ncol(feature.set), "\n",
        file = paste0(res.dir, "/log.txt"), append = TRUE)
  } else if(featset.name == "union"){
    # select union descriptor
    union.feat <- get.union.feat(final.comp, prot)
    feature.set <- descriptors[, union.feat]
    cat("  Feature set name:", featset.name, ", number:", ncol(feature.set), "\n",
        file = paste0(res.dir, "/log.txt"), append = TRUE)
  } else { stop("Selected feature set name is not available") }
  
  # Load models from RData file
  load(file.path(MODELS.DIR, prot, paste0(featset.name, "_models_", prot, ".RData")))
  
  # Predict activity with best model
  switch(bmodel,
         lasso={ current.model <- lasso.model
         pred <- predict(lasso.model, newx = as.matrix(feature.set), type = "class")
         pred <- factor(ifelse(pred == "1", "ACTIVE", "INACTIVE"), levels = c("ACTIVE", "INACTIVE"))
         },
         svm={ current.model <- svm.model
         pred <- predict(current.model, feature.set)
         },
         decision_tree={ current.model <- c50.model
         pred <- predict(current.model, feature.set)
         },
         xgboost={ current.model <- xgboost.model
         pred <- predict(xgboost.model, data.matrix(feature.set))
         pred <- factor(ifelse(pred > 0.5, "ACTIVE", "INACTIVE"), levels = c("ACTIVE", "INACTIVE"))
         },
         logistic={ current.model <- logistic.model
         pred <- predict(logistic.model, feature.set, type="response")
         pred <- factor(ifelse(pred > 0.5, "INACTIVE", "ACTIVE"), levels = c("ACTIVE", "INACTIVE"))
         },
         random_forest={ current.model <- rf.model
         pred <- predict(current.model, feature.set)
         },
         neuralnet={ current.model <- neuralnet.model
         pred <- predict(neuralnet.model, feature.set)
         pred <- factor(ifelse(pred[,1] > 0.5, "ACTIVE", "INACTIVE"), levels = c("ACTIVE", "INACTIVE"))
         },
         knn={ current.model <- knn.model
         pred <- predict(current.model, feature.set)
         },
         elasticnet={ current.model <- elnet.model
         pred <- predict(current.model, feature.set)
         },
         naive={ current.model <- naive.model
         pred <- predict(current.model, feature.set)
         },
         cart={ current.model <- cart.model
         pred <- predict(cart.model, feature.set, type = "class")
         },
         ridge={ current.model <- ridge.model
         pred <- predict(ridge.model, newx = as.matrix(feature.set), type = "class")
         pred <- factor(ifelse(pred == "1", "ACTIVE", "INACTIVE"), levels = c("ACTIVE", "INACTIVE"))
         },
         { current.model <- NULL }
  )
  
  # Remove unused models for next execution
  rm(current.model, naive.model, logistic.model, svm.model, c50.model,
     rf.model, neuralnet.model, xgboost.model, knn.model,
     cart.model, lasso.model, ridge.model, elnet.model)
  
  return(data.frame(class=pred, Name=rownames(descriptors),prot))
})

# Collect and save prediction results
kuala.results <- do.call(what = rbind, args = res)
write.table(kuala.results, file.path(res.dir,"kuala_predicted_ligands_activity.txt"),
            sep="\t", quote=FALSE, row.names = FALSE)
