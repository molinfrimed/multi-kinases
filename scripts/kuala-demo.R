options(stringsAsFactors = FALSE)

# Load functions
source("Script/model_functions.R")

raw.descrs <- read.delim("Version_2/feat_selec_003_10.csv")
descrs.all <- cbind(raw.descrs[,c(1:3)], scale(raw.descrs[,-c(1:3)]))


# Compare ChEMBL ids
pkidb <- read.delim("Data_DB/Data_DB.tsv")
pkidb.saved <- pkidb[which(!pkidb$CHEMBL_ID %in% descrs.all$Chembl_ID), ]
# write.table(pkidb.saved, file = "Data_DB/new_ligands.txt", quote = FALSE, sep = "\t", row.names = FALSE)


## exclude target id  == nan
pkidb.ligands <- pkidb.saved[which(!pkidb.saved$targets_id == "nan"), "Name"]

# Load computed descriptors
descriptors <- read.delim("ligands_descriptors.tsv")
# raw.descriptors <- raw.descriptors[which(raw.descriptors$Name %in% pkidb.ligands),]

# res <- apply(descriptors[,-1], 2, function(feat){
#   na.counts <- sum(is.na(feat)) # Count NA values
#   all.zeros <- all(feat == 0) # all zeros
#   zero.var <- sd(feat, na.rm = TRUE) # zero-variance
#   return(data.frame(na.counts, all.zeros, zero.var))
# })
# stats <- do.call(rbind, res)
# na.feat <- rownames(stats)[which(stats$na.counts > 0)]
# rm(stats)
# which(apply(descriptors, 2, is.infinite))

# descriptors <- data.frame(scale(raw.descriptors[,-1]))
# rownames(descriptors) <- raw.descriptors$Name
# write.table(descriptors, "github/ligands_descriptors.tsv", quote = FALSE, sep = "\t")

# replace NaN values with 0
descriptors[,names(which(apply(descriptors, 2, function(col){
  all(is.nan(col))
})))] <- 0
# new.descriptors <- merge(raw.descriptors, pkidb.saved, by="Name")

##############

# Load best models and set metric
best.models <- read.delim("best_mfs_kuala.tsv")
metric <- "F1"
res.dir <- "kuala_results"

# Create complete set of ligands
# ligands <- unique(descrs.all$CHEMBL_ID)
proteins <- unique(best.models$prot)

# Load feature set comparison
final.comp <- read.delim(file = "feat_sel_df.csv")

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
  load(paste0("Version_2/RESULTS/", prot, "/", featset.name, "_models_", prot, ".RData"))
  
  # Select best model
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
  
  # Remove models
  rm(current.model, naive.model, logistic.model, svm.model, c50.model,
     rf.model, neuralnet.model, xgboost.model, knn.model,
     cart.model, lasso.model, ridge.model, elnet.model)
  
  return(data.frame(class=pred, Name=rownames(descriptors),prot))
})

pkidb.results <- do.call(rbind, res)
write.table(pkidb.results, "Version_2/validation_results/pkidb_results.tsv", sep="\t", quote=FALSE,
            row.names = FALSE)

convert.df <- do.call(rbind, lapply(1:nrow(pkidb.saved), function(i){
  data.frame(Name=pkidb.saved$Name[i],
             uniprot=unlist(strsplit(pkidb.saved$targets_id[i], "|", fixed = TRUE)))
}))

classif.annots <- merge(pkidb.results, convert.df, by="Name")
validation.results <- classif.annots[classif.annots$prot == classif.annots$uniprot, ]
write.table(validation.results, "Version_2/validation_results/pkidb_results_match.tsv", sep="\t", quote=FALSE,
            row.names = FALSE)




