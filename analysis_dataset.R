#!/usr/bin/env Rscript

library("argparse")
suppressPackageStartupMessages(library("argparse"))

help_message <-
"usage: Rscript analysis_dataset.R
	<--cohort A>
	<--dspath /home/kenneth/out/UKBB/covid19/cohorts.RData>
	<--metapath /home/kenneth/data/UKBB/covid19/meta/5Jan2021>
	<--pval FALSE>
	<--cpu 40>
	<--outpath /home/kenneth/out/UKBB/covid19/run_03Nov2020>
\nRequired:
	--cohort    Cohort (A,B,C,D) to be analysed
	--dspath    RData file path storing cohorts A-D
	--metapath  Path to the metadata about the selected traits
	--pval      Enable the p-value estimation
	--cpu       Total number of Cores used
	--outpath   Output folder
	--help | -h Display this help message\n"

parser <- ArgumentParser()
parser$add_argument("--cohort", action="store", type="character", default="A", dest="cohort")
parser$add_argument("--dspath", action="store", dest="dspath")
parser$add_argument("--metapath", action="store", dest="metapath", default=NULL)
parser$add_argument("--pval", action="store", type="logical", default=FALSE, dest="enable_pval")
parser$add_argument("--cpu", action="store", type="integer", dest="cpu", default=40)
parser$add_argument("--outpath", action="store", dest="outpath")

capture <- commandArgs(trailingOnly = TRUE)
help <- (sum(c("--help", "-h") %in% capture) >= 1)
if (help) {
	cat(help_message)
	quit()
}


# Parse user arguments
args <- parser$parse_args()

total_cores=args$cpu

# Load the cohort A,B,C,D from RData files
load(args$dspath)

dim(df.a)
dim(df.b)
dim(df.c)
dim(df.d)
# ________________ One-hot encoding for the cohort________________________

library(pROC)
library(plyr)
library(dplyr)
library(tidyr)
library(pscl)
library(xgboost)
library(caret)
library(e1071)
library(auctestr)	# For calculation of SE of AUC

# Model Building (XGBoost for all Traits (CAT + QTL))

# Fold 1-7 seed = c(7, 12, 230, 103, 88, 38, 32)
set.seed(32)

# *******************************************************
# * Change the Cohort to A or B or C                    *
# * Then re-run the block of code below this line       *
# *******************************************************
cohort = args$cohort

# Build model based on specified cohort
if(cohort == "A") {
  df.final <- df.a
} else if(cohort == "B") {
  df.final <- df.b
} else if(cohort == "C") {
  df.final <- df.c
} else if(cohort == "D") {
  df.final <- df.d
}

# Read the list of categorical and quantitative features to be included in the model buildling
if(!is.null(args$metapath)) {
	selected_cat_traits <- scan(paste0(args$metapath, "/selected_cat_features.txt"), character(), sep="\n")
	selected_qtl_traits <- scan(paste0(args$metapath, "/selected_qtl_features.txt"), character(), sep="\n")
}
print(paste0("CAT traits (", length(selected_cat_traits), "), QTL traits (", length(selected_qtl_traits), ")"))


#______________________________________________________________________
#      this line changed to avoid collinear variables (9 Sep 2020)
#______________________________________________________________________
# One Hot Enocde the Categorical predictor
dmy <- dummyVars(" ~ .", data = df.final[selected_cat_traits], sep=".", fullRank=TRUE) 
#https://amunategui.github.io/dummyVar-Walkthrough/)
df.cat_traits.hotenc <- data.frame(predict(dmy, newdata = df.final[selected_cat_traits]))

# Replace the categorical features by the Hot encoding features
df.final <- cbind(df.cat_traits.hotenc, df.final[,c(selected_qtl_traits, "outcome")])

dim(df.final)
dim(subset(df.final, outcome=="case"))
dim(subset(df.final, outcome=="ctrl"))


# Randomly shuffle the data
# https://gist.github.com/duttashi/a51c71acb7388c535e30b57854598e77  

N_folds <- 5
df.final.shuf <- df.final[sample(nrow(df.final)),]
folds <- cut(seq(1,nrow(df.final.shuf)),breaks=N_folds,labels=FALSE) #Create N equally size folds

#Perform k-fold cross validation
testData_list <- vector(mode = "list", length = N_folds)
trainData_list <- vector(mode = "list", length = N_folds)
for (i in 1:N_folds) {

  #Segement your data by fold using the which() function
  testIndexes <- which(folds==i,arr.ind=TRUE)
  testData_list[[i]] <- df.final.shuf[testIndexes, ]
  trainData_list[[i]] <- df.final.shuf[-c(testIndexes), ]
  
  #Use the test and train data partitions however you desire...
}


# Tune the hyperparameters of XGBoost model for each training set
params_best <- list()
for(i in 1:N_folds){
  df.train = trainData_list[[i]]
  
  # Set label
  xgb.train.label <- as.integer(df.train$outcome == "case")
  
  # Convert dataframe to numeric type
  df.train.num <- sapply(df.train[,1:ncol(df.train)-1], as.numeric)
  
  # Convert dataframe to matrix type
  dtrain <- xgb.DMatrix(data=as.matrix(df.train.num), label=xgb.train.label)
  
  # HyperParameters tuning by formal grid search
  # https://stats.stackexchange.com/questions/171043/how-to-tune-hyperparameters-of-xgboost-trees
  # XGBoost task parameters
  nrounds <- 80
  folds <- 5
  obj <- 'binary:logistic'
  
  library("NMOF")
  # Parameter grid to search (please change as needed) 
  # Cohort A
  paramsA <- list(
    eta = c(0.15,0.2,0.25),
    max_depth = 2:4,
    # max_delta_step = c(0,1),
    gamma= c(0.5,1,1.5),
    lambda= c(3,5,10), #L2 penalty
    alpha= c(8,10,15), #L1 penalty
    nthread=2,
    subsample = 1
  )
  
  
  # Cohort B
  ## Logloss
  paramsB <- list(
    eta = c(0.05, 0.1,0.2),
    max_depth = 2:3,
    # max_delta_step = c(0,1),
    gamma= c(0.1, 1, 5),
    lambda= c(1,5,10), #L2 penalty
    alpha= c(6,12,15), #L1 penalty
    nthread=2,
    subsample = 1
  )
  
  # Cohort C
  ctrl_case_ratio <- length(which(df.c$outcome == "ctrl"))/length(which(df.a$outcome == "case"))
  paramsC <- list(
    eta = c(0.2),
    tree_method = c("hist"),
    max_depth = c(5,6),
    max_delta_step = c(0),
    min_child_weight = c(5,8),
    gamma= c(0.1, 0.01),
    lambda= c(5), #L2 penalty
    alpha= c(15), #L1 penalty
    scale_pos_weight = sqrt(ctrl_case_ratio),
    nthread=4,
    subsample = 1
  )

  ctrl_case_ratio <- length(which(df.c$outcome == "ctrl"))/length(which(df.b$outcome == "case"))
  # Cohort D
  paramsD <- list(
    eta = c(0.2),
    tree_method = c("hist"),
    max_depth = c(5,6),
    max_delta_step = c(0),
    min_child_weight = c(5,8),
    gamma= c(0.1, 0.01),
    lambda= c(5), #L2 penalty
    alpha= c(15), #L1 penalty
    scale_pos_weight = sqrt(ctrl_case_ratio),
    nthread=4,
    subsample = 1
  )
    
  if(cohort == "A"){
    params <- paramsA
  } else if (cohort == "B"){
    params <- paramsB
  } else if (cohort == "C"){
    params <- paramsC
  } else if (cohort == "D"){
    params <- paramsD
  }
  
  library(parallel) 
  # Table to track performance from each worker node
  res <- data.frame()
  
  # Simple cross validated XGboost training function (returning minimum error for grid search)
  # Here we only use the training set for finding best tuning parameters
  xgbCV <- function (search_params) {
    fit <- xgb.cv(
  	  data = as.matrix(df.train.num),
  	  label = xgb.train.label,
      param = search_params, 
      missing = NA, 
      nfold = folds, 
  	  objective = obj,
      prediction = FALSE,
      early_stopping_rounds = 20,
      maximize = TRUE,
  	  metrics = c("auc"), #c("logloss"),
      nrounds = nrounds,
  	  stratified=T, 
  	  print_every_n=10
    )
    
    fit_df = fit$evaluation_log
    fit_params = fit$params
    rounds <- nrow(fit_df)
    metric <- ifelse(fit_params$eval_metric=="auc","test_auc_mean", "test_logloss_mean")
    # metric = "test_auc_mean"  #"test_logloss_mean"    #paste('test.',eval,'.mean',sep='')
    # metric = "test_logloss_mean"
    if (metric =="test_logloss_mean") {
      idx <- which.min(fit_df[,fit_df[[metric]]])
      retval <- fit_df[idx,][[metric]]
    }
    if (metric =="test_auc_mean") {
      idx <- which.max(fit_df[,fit_df[[metric]]])
      retval <- 1.0 - fit_df[idx,][[metric]]
    }
    val <- fit_df[idx,][[metric]]
    
    res <<- rbind(res, c(idx,val,metric,rounds))

    colnames(res) <<- c('idx','val','metric','rounds')
    return(retval)
    
  }
  
  # Find minimal testing error in parallel
  cl <- makeCluster(round(total_cores/2)) 
  
  library(doParallel) 
  registerDoParallel(cl)
  clusterExport(cl, c("xgb.cv","df.train.num", "xgb.train.label", 
                      "nrounds", "obj","res","eval","folds","params"))
  
  t1= proc.time()
  sol <- gridSearch(
    fun = xgbCV,
    levels = params,
    method =  'snow', #"multicore",
    cl = cl,
    #mc.control = list(mc.cores=8), 
    keepNames = TRUE,
    asList = TRUE
  )
  proc.time()-t1
  
  # Combine all model results
  comb=clusterEvalQ(cl,res)
  results <- ldply(comb,data.frame)
  stopCluster(cl)
  params_best <- append(params_best, list(sol$minlevels))
  # results
}
# ____________________end of parameter tuning_________________________


# ______________Build the XGB Model for K-fold CV_____________________
df.out.final <- data.frame()
pred.auc <- c()
xgb.models <- list()

for(i in 1:N_folds){
  df.train = trainData_list[[i]]
  df.test = testData_list[[i]]
  
  # Further split train set to train and validation set
  train1.idx <- sample(nrow(df.train), 4/5 * nrow(df.train))
  df.train1 <- df.train[train1.idx,]
  df.train2 <- df.train[-train1.idx,]

  # Set label
  xgb.train.label <- as.integer(df.train$outcome == "case")
  xgb.train1.label <- as.integer(df.train1$outcome == "case")
  xgb.train2.label <- as.integer(df.train2$outcome == "case")
  xgb.test.label <- as.integer(df.test$outcome == "case")

  # Convert dataframe to numeric type
  df.train.num <- sapply(df.train[,1:ncol(df.train)-1], as.numeric)
  df.train1.num <- sapply(df.train1[,1:ncol(df.train1)-1], as.numeric)
  df.train2.num <- sapply(df.train2[,1:ncol(df.train2)-1], as.numeric)
  df.test.num <- sapply(df.test[,1:ncol(df.test)-1], as.numeric)

  # Convert dataframe to matrix type
  dtrain <- xgb.DMatrix(data=as.matrix(df.train.num), label=xgb.train.label)
  dtrain1 <- xgb.DMatrix(data=as.matrix(df.train1.num), label=xgb.train1.label)
  dtrain2 <- xgb.DMatrix(data=as.matrix(df.train2.num), label=xgb.train2.label)
  dtest <- xgb.DMatrix(data=as.matrix(df.test.num), label=xgb.test.label)
  watchlist <- list(train=dtrain1, test=dtrain2)
  
  # Find the optimium iteration
  xgbcv <- xgb.cv(params=params_best[[i]], data=dtrain1, label = NULL,
                  nrounds=nrounds, nfold=5, showsd=T, 
                  stratified=T, print,every.n=10, 
                  early_stopping_rounds=20, maximize=T, 
                  eval_metric="auc")
  
  # Build model
  xgb <- xgb.train(data=dtrain1, watchlist=watchlist, 
                   params=params_best[[i]],# nthread=20, 
                   nrounds=xgbcv$best_iteration,
                   objective="binary:logistic", 
                   eval_metric="auc", verbose=2)
  xgb.models <- append(xgb.models, list(xgb))
  
  # Evaluate testing AUC
  test.xgb <- predict(xgb, as.matrix(df.test.num))
  
  pred.roc <- roc(xgb.test.label, test.xgb, direction="<",
                  plot=TRUE, percent=TRUE, ci=TRUE, 
                  ci.alpha=0.95, grid=TRUE, 
                  print.auc=TRUE, show.thres=TRUE)
  pred.auc <- c(pred.auc, auc(pred.roc))

  # Create dataframe for subjects in test set and sorted by predicted clinical risk
  df.out.test <- df.test["outcome"]
  df.out.test$clinical_risk <- test.xgb
  df.out.test$eid <- row.names(df.out.test)
  
  # Append the subjects with clinical risk to dataframe
  df.out.final <- rbind(df.out.final, df.out.test)

}

# Print the Average AUC performance metric
print("The 5 fold AUC:\n")
print(pred.auc)
print(paste0("Avg Testing AUC of 5 fold CV set is ", round(mean(pred.auc),1), "%"))

# Calculate the SE and 95% CI of AUC
auc.summary <- data.frame(matrix(ncol = 5, nrow = 0)) # Fold, AUC, SE
for (i in 1:N_folds) {

  # Calculate the number of case and control in testing dataset
  df <- testData_list[[i]]
  n_total <- nrow(df)
  n_case <- length(which(df$outcome == "case"))
  n_ctrl <- nrow(df) - n_case
      
  # Calculate the Standard Error of AUC
  auc <- pred.auc[i]/100
  auc.se <- se_auc(auc, n_case, n_ctrl)
        
  # Evaluate the mean error for 95% CI of the AUC
  me <- qnorm(0.975) * auc.se
  auc.lower_ci <- auc - me
  auc.upper_ci <- auc + me
	    
  print(paste0("Fold ", i, ": AUC=", auc, ", SE=", auc.se, 
               ", 95% CI=[", auc.lower_ci, ",", auc.upper_ci, "]"))
  
  auc.summary <- rbind(auc.summary,
		       list(i, auc, auc.se, auc.lower_ci, auc.upper_ci))
}
names(auc.summary) <- c("fold", "AUC", "AUC.SE", "AUC.UPPER_CI", "AUC.LOWER_CI")

# Estimate the Average AUC SE and 95% CI
auc.avg <- mean(auc.summary$AUC)
auc.avg.se <- sqrt(sum(auc.summary$AUC.SE^2)/(N_folds^2))
me <- qnorm(0.975) * auc.avg.se
auc.avg.lower_ci <- auc.avg - me
auc.avg.upper_ci <- auc.avg + me

auc.summary <- rbind(auc.summary,
		     list("Avg", auc.avg, auc.avg.se,
			  auc.avg.lower_ci, auc.avg.upper_ci))

# Save the AUC statistical summary to file
auc.summary[,-1] <- round(auc.summary[,-1], 4)
write.table(auc.summary, row.names=FALSE, quote=FALSE, sep="\t",
            file=paste0(args$outpath, "/cohort", cohort, "_auc_summary.tsv"))


# Save the subjects with sorted clinical risk to file
df.out.final.sorted <- df.out.final[order(-df.out.final[,2]),]
write.table(df.out.final.sorted[,c("eid","clinical_risk","outcome")],
	    paste0(args$outpath, "/cohort", cohort, "_clinical_risk.tsv"),
            row.names=FALSE, quote=FALSE, sep="\t")

save(trainData_list, testData_list, xgb.models, file=paste0(args$outpath, "/xgb_model/cohort", cohort, "_XGBModels.RData"))


# Save the xgb_model and test_auc to file
setwd(paste0(args$outpath, "/xgb_model"))
for (i in 1:N_folds) {
  best_xgb_model <- xgb.models[[i]]
  test.auc <- pred.auc[i]
  
  # Save the best model
  save(best_xgb_model, test.auc, file=paste0("xgb_model_cohort_", cohort, "_fold_", i, ".Rdata"))
}


#________________________________________________________________
#** now we compute average shapley value over all folds ** 
#________________________________________________________________
library(SHAPforxgboost)
setwd(paste0(args$outpath, "/xgb_model"))

dataX <- as.matrix(rbind(df.train.num, df.test.num))
rownames(dataX) <- rownames(rbind(df.train, df.test))
save(dataX, file=paste0("dataX_trainTestcomb_cohort", cohort, ".Rdata"))

list.shap_val = list()

#dataX is fixed 
for (i in 1:N_folds) {
  xgb_model <- xgb.models[[i]]
  list.shap_val[[i]] <- predict(xgb_model, dataX, predcontrib = TRUE, approxcontrib = F)
}

sum.shap_val = 0 
for (i in 1:N_folds) {
  sum.shap_val <- list.shap_val[[i]]  + sum.shap_val
}

# Save the mean shapley value
mean.shap_val <- sum.shap_val/N_folds
save(mean.shap_val, file=paste0("Mean_shap_cohort", cohort,".Rdata"))


#** sort the variables by mean |shap| **
ind_bias = which(colnames(mean.shap_val) == "BIAS")
mat_meanAbsShap = data.frame( Feature = colnames(mean.shap_val[,-ind_bias]) , 
			                 meanAbsShap = colMeans( abs(mean.shap_val[,-ind_bias]) ) 
					            )

mat_meanAbsShap$Feature = as.character(mat_meanAbsShap$Feature)
mat_meanAbsShap <- arrange(mat_meanAbsShap, -meanAbsShap)


# ________________________ p value calculation________________________

# Compute the p-values for the varImp (as quantified by mean abs shapley value)
if (args$enable_pval) {

  set.seed(32)
  nrounds = 80
  no_perm = 500  # For cohort C and D, please change to 100 to save time
  best_fold_idx = which.max(pred.auc)
  best_params <- params_best[[best_fold_idx]]
  best_params$nthread = 4

  # Extract the feature values only
  df.features <- subset(df.final, select=-c(outcome))

  # Build the XGB Model for no_perm times
  df.meanAbsShap.null <- data.frame(matrix(ncol = ncol(df.features), nrow = 0))
  for(i in 1:no_perm){
    outcome <- sample(df.final$outcome, replace=F)
    df.train <- cbind(df.features, outcome)
    dataX <- as.matrix(df.train)

    # Further split train set to train and validation set
    train1.idx <- sample(nrow(df.train), 4/5 * nrow(df.train))
    df.train1 <- df.train[train1.idx,]
    df.train2 <- df.train[-train1.idx,]

    # Set label
    xgb.train.label <- as.integer(df.train$outcome == "case")
    xgb.train1.label <- as.integer(df.train1$outcome == "case")
    xgb.train2.label <- as.integer(df.train2$outcome == "case")

    # Convert dataframe to numeric type
    df.train.num <- sapply(df.train[,1:ncol(df.train)-1], as.numeric)
    df.train1.num <- sapply(df.train1[,1:ncol(df.train1)-1], as.numeric)
    df.train2.num <- sapply(df.train2[,1:ncol(df.train2)-1], as.numeric)

    # Convert dataframe to matrix type
    dtrain <- xgb.DMatrix(data=as.matrix(df.train.num), label=xgb.train.label)
    dtrain1 <- xgb.DMatrix(data=as.matrix(df.train1.num), label=xgb.train1.label)
    dtrain2 <- xgb.DMatrix(data=as.matrix(df.train2.num), label=xgb.train2.label)
    watchlist <- list(train=dtrain1, test=dtrain2)

    # Find the best number of iterations
    xgbcv <- xgb.cv(params=best_params, data=dtrain1, label = NULL,
                  nrounds=nrounds, nfold=5, showsd=T, 
                  stratified=T, print,every.n=20, 
                  early_stopping_rounds=20, maximize=T, 
                  eval_metric="auc", verbose=0)
  
    # Build model
    xgb_model_shuffled_outcome <- xgb.train(data=dtrain1, watchlist=watchlist, 
                   params=best_params, # nthread=20, 
                   nrounds=xgbcv$best_iteration,
                   objective="binary:logistic", 
                   eval_metric="auc", verbose=0)

    # Predict the shapley value
    mean.shap_val.null <- predict(xgb_model_shuffled_outcome, dtrain, predcontrib=TRUE, approxcontrib=F)  
  
    # Append the colMeans of meanAbsShap of all features to dataframe
    ind_bias <- which(colnames(mean.shap_val.null) == "BIAS")
    df.meanAbsShap.null <- rbind(df.meanAbsShap.null, colMeans(abs(mean.shap_val.null[,-ind_bias])))
  }
  colnames(df.meanAbsShap.null) <- as.character(colnames(mean.shap_val.null[,-ind_bias]))

  #match the feature names in the same order 
  id_mat = match(mat_meanAbsShap$Feature,  table=colnames(df.meanAbsShap.null))
  df.meanAbsShap.null2 <- df.meanAbsShap.null[,id_mat]

  pval = c()
  no_perm = nrow(df.meanAbsShap.null2)
  for (i in 1:ncol(df.meanAbsShap.null2)) {
    pval[i] <-(sum(df.meanAbsShap.null2[,i]>=mat_meanAbsShap$meanAbsShap[i] ) + 1) / (no_perm+1)
  }

  df.meanAbsShap.obs.out <- data.frame(Feature=mat_meanAbsShap$Feature,
                                     meanAbsShap=mat_meanAbsShap$meanAbsShap,
                                     pval=pval)

  setwd(paste0(args$outpath, "/pvalue"))
  # Save the Null and Observed distribution of mean absolute Shapley Value
  save(no_perm, df.meanAbsShap.null2, mat_meanAbsShap,
       file=paste0("Mean_shap_cohort", cohort, "_nullvstest.Rdata"))

  write.table(t(df.meanAbsShap.obs.out), col.names=FALSE, row.names=FALSE, quote=FALSE, sep=",",
       file=paste0("Mean_shap_cohort", cohort, "_pval.csv"))

}

