#!/usr/bin/env Rscript

library("argparse")
suppressPackageStartupMessages(library("argparse"))

help_message <-
  "usage: Rscript analysis_dataset.R
        <--cohort A>
        <--dspath /home/kenneth/out/UKBB/covid19/run_03Nov2020/impute_with_replace/xgb_model>
	<--naming /home/kenneth/data/UKBB/covid19/meta/names_lookup.csv>
        <--outpath /home/kenneth/out/UKBB/covid19/run_03Nov2020/impute_with_replace>
\nRequired:
        --cohort    Cohort (A,B,C,D) to be analysed
        --dspath    Folder storing training data and XGboost built model
	--naming    A dictionary file (covariates names => label in graphs)
        --outpath   Output folder
        --help | -h Display this help message\n"

parser <- ArgumentParser()
parser$add_argument("--cohort", action="store", type="character", default="A", dest="cohort")
parser$add_argument("--dspath", action="store")
parser$add_argument("--naming", action="store")
parser$add_argument("--outpath", action="store")

capture <- commandArgs(trailingOnly = TRUE)
help <- (sum(c("--help", "-h") %in% capture) >= 1)
if (help) {
  cat(help_message)
  quit()
}

# Parse user arguments
args <- parser$parse_args()


#___________________________________________________________________________________
# changed 10 Sep 2020 
# in view of the instability of shap and varImp, we take average over 5-folds 
#_____________________________________________________________________________________
#devtools::install_github("liuyanguu/SHAPforxgboost")
library(SHAPforxgboost)
library(xgboost)
library(dplyr)
library(data.table)

#assumed to be followed by fold number and .Rdata
N_folds = 5 
cohort = args$cohort

dataX_path = paste0(args$dspath, "/dataX_trainTestcomb_cohort", cohort, ".Rdata")
xgbmodel_path = paste0(args$dspath, "/xgb_model_cohort_", cohort, "_fold_")
out_path = args$outpath # "/home/kenneth/out/UKBB/covid19/run_03Nov2020/impute_with_replace"

df.nlookup <- read.csv(args$naming, header = TRUE, stringsAsFactors=FALSE, quote="\"")

# Set the output path of saved graph and Rdata
setwd(paste0(out_path, "/xgb_model"))

# Load the training data "dataX"
load(dataX_path)

# Variables initialization
list.varimp <- list()
no_features = ncol(dataX)
feature_names = colnames(dataX)
xgb.models = list()

# Calculate the average shapley value for each features
for (i in 1:N_folds) {
  load(paste0(xgbmodel_path, i, ".Rdata"))
  xgb.models[[i]]  <-  best_xgb_model
  
  # Save the best model
  list.varimp[[i]] <- xgb.importance(feature_names= feature_names , model=best_xgb_model)
  
  # features that are not matched are assumed to have importance of zero
  match_ind = match(list.varimp[[i]]$Feature, table=feature_names )
  not_matched = setdiff(1:no_features, match_ind)
  zero_imp_mat = matrix(0, nrow=length(not_matched), ncol=3)
  zero_mat = data.frame(feature_names[not_matched], zero_imp_mat)
  colnames(zero_mat) <- c("Feature","Gain", "Cover", "Frequency")
  list.varimp[[i]]  = rbind(list.varimp[[i]], zero_mat)
  list.varimp[[i]] = data.frame(list.varimp[[i]])
  
  # sort by alphabetical order of Feature column so all matrices match
  list.varimp[[i]]$Feature <- as.character(list.varimp[[i]]$Feature)
  list.varimp[[i]] = list.varimp[[i]] [order(list.varimp[[i]]$Feature),]
}

sum.varimp = 0
for (i in 1:N_folds) {
  sum.varimp <- list.varimp[[i]][,2:4] + sum.varimp
}

mean_list.varimp = data.frame( Feature=list.varimp[[1]]$Feature, sum.varimp/N_folds)
mean_list.varimp  = arrange(mean_list.varimp, -Gain)
mean_list.varimp = as.data.table(mean_list.varimp)

# Save the varImp to file
save(mean_list.varimp, file=paste0("mean_list.varimp_cohort", cohort,".Rdata"))

# Replace the covariates name to full description for plotting graphs
mean_list.varimp$Feature = df.nlookup$desc[match(mean_list.varimp$Feature, df.nlookup$predictor)]

#_______________________________________________
#**Mean XGboost VarImp plot**
#______________________________________________
library(ggplot2)


# Plot the varImp (top 30 features)
p <- xgb.ggplot.importance(importance_matrix=mean_list.varimp[1:30])
p+theme_bw()+
  theme(axis.text=element_text(face="plain", angle=30, size="11", color="black"),
	panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
	legend.position="none")+
  labs(x="Risk Factor", y="XGBoost VarImp", title="")

ggsave(device="png", 
       file=paste0("Mean_xgb.importance_cohort", cohort, "_varimp_plot_top30", ".png"),
       width=6, height=7, unit="in", dpi=300)
dev.off()

# Plot the varImp (top 15 features)
p <- xgb.ggplot.importance(importance_matrix=mean_list.varimp[1:15])
p+theme_bw()+
  theme(axis.text=element_text(face="plain", angle=30, size="11", color="black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position="none")+
  labs(x="Risk Factor", y="XGBoost VarImp", title="")

ggsave(device="png",
       file=paste0("Mean_xgb.importance_cohort", cohort, "_varimp_plot_top15", ".png"),
       width=6, height=7, unit="in", dpi=300)
dev.off()


#________________________________________________________________
#** now we compute average shapley value over all folds ** 
#________________________________________________________________


setwd(paste0(out_path, "/shapley"))

#dataX is fixed 
list.shap_val = list()
for (i in 1:N_folds) {
  list.shap_val[[i]] <- predict(xgb.models[[i]], dataX, predcontrib = TRUE, approxcontrib = F)
}

sum.shap_val = 0 
for (i in 1:N_folds) {
  sum.shap_val <- list.shap_val[[i]]  + sum.shap_val
}

mean.shap_val <- sum.shap_val/N_folds
save(mean.shap_val, file=paste0("mean.shap_val_cohort", cohort,".Rdata"))

# Replace the covariate name by full description
ind_bias = which(colnames(mean.shap_val) == "BIAS")
colnames(mean.shap_val)[-ind_bias] <- df.nlookup$desc[match(colnames(mean.shap_val[,-ind_bias]), df.nlookup$predictor)]

dataX.rename <- dataX
colnames(dataX.rename) <- df.nlookup$desc[match(colnames(dataX.rename), df.nlookup$predictor)]

#** sort the variables by mean |shap| **
mat_meanAbsShap = data.frame( Feature = colnames(mean.shap_val[,-ind_bias]) , 
                              meanAbsShap = colMeans( abs(mean.shap_val[,-ind_bias]) ) 
)
mat_meanAbsShap$Feature = as.character(mat_meanAbsShap$Feature)
mat_meanAbsShap <- arrange(mat_meanAbsShap, -meanAbsShap) 

#_______________________________________________
#**Mean Shapley VarImp plot**
#______________________________________________

setwd(paste0(out_path, "/shapley"))

# SHAP VarImp Plot
pdf(paste0("Mean_shap_cohort", cohort, "_varimp_plot_top30.pdf"),
    width=7, height=9.5, paper="a4")
xgb.plot.shap(data=dataX.rename,
              shap_contrib=mean.shap_val,
              features=mat_meanAbsShap$Feature[1:30],
              n_col=3)
dev.off()

pdf(paste0("Mean_shap_cohort", cohort, "_varimp_plot_top15.pdf"),
    width=7, height=4.75, paper="a4")
xgb.plot.shap(data=dataX.rename,
              shap_contrib=mean.shap_val,
              features=mat_meanAbsShap$Feature[1:15],
              n_col=3)
dev.off()

pdf(paste0("Mean_shap_cohort", cohort, "_varimp_plot_top6.pdf"),
    width=7, height=2.85, paper="a4")
xgb.plot.shap(data=dataX.rename,
              shap_contrib=mean.shap_val,
              features=mat_meanAbsShap$Feature[1:6],
              n_col=3)
dev.off()

# SHAP VarImp Summary Plot
if (cohort == "C" | cohort =="D") {
  subsamp_N = 10000
  ind_subsamp = sample(nrow(dataX), size=subsamp_N)
  shap.score = data.frame(mean.shap_val[,-ind_bias][ind_subsamp,])
  dataX.1 <- dataX.rename[ind_subsamp,]
} else {
  shap.score = data.frame(mean.shap_val[,-ind_bias])
  dataX.1 <- dataX.rename
}
pdf(paste0("Mean_shap_Overall_cohort", cohort, "_varimp_plot_top30.pdf"),
    width=7, height=9.5, paper="a4")
shap.plot.summary.wrap2(shap_score= shap.score, 
                        X= dataX.1, 
                        top_n=30,
                        dilute = FALSE)
dev.off()

pdf(paste0("Mean_shap_Overall_cohort", cohort, "_varimp_plot_top15.pdf"),
    width=7, height=4.75, paper="a4")
shap.plot.summary.wrap2(shap_score= shap.score, 
                        X= dataX.1, 
                        top_n=15,
                        dilute = FALSE)
dev.off()


#_______________________________________________
#** Shapley interaction plot**
#______________________________________________

shap_int.absmean = list()

#if too time-consuming, may randomly select some individuals
dataX.1 <- dataX
ind_subsamp <- NULL
if (cohort == "C" | cohort == "D") {
  subsamp_N = 50000
  ind_subsamp = sample(nrow(dataX), size=subsamp_N)
  dataX.1 <- dataX[ind_subsamp,]
}


for (i in 1:N_folds) {
  shap_int <- predict(xgb.models[[i]], 
                      dataX.1, 
                      predinteraction = TRUE)
  
  abs_shap_int <- abs(shap_int)
  
  #taking the mean over dim 2 and 3, ie mean of each subject (dim 1 is the subject)
  shap_int.absmean[[i]] <- apply(abs_shap_int , c(2,3), sum)
  
  #remove bias row and column
  shap_int.absmean[[i]] <- shap_int.absmean[[i]] [-nrow(shap_int.absmean[[i]]), -ncol(shap_int.absmean[[i]])]
  
}

sum.int_val = 0 
for (i in 1:N_folds) {
  sum.int_val <- shap_int.absmean[[i]]  + sum.int_val
}

mean.int_val = sum.int_val/N_folds

library(dplyr)
library(reshape2)
mean.int_val[upper.tri(mean.int_val)] <- NA
mean.int_val_list <- melt(mean.int_val, na.rm=TRUE)
ind_interact = which(mean.int_val_list$Var1 != mean.int_val_list$Var2)
mean.int_val_list_diffTraits = arrange(mean.int_val_list[ind_interact,], -value) 

setwd(paste0(out_path, "/shapley"))
save(ind_subsamp, mean.int_val, mean.int_val_list_diffTraits, file= paste0("Interaction_Shap", cohort, ".Rdata"))


mean.int_val_list_diffTraits$Var1 = as.character(mean.int_val_list_diffTraits$Var1)
mean.int_val_list_diffTraits$Var2 = as.character(mean.int_val_list_diffTraits$Var2)


# Replace the covarite name by full description
mean.int_val_list_diffTraits$Var1 <- df.nlookup$desc[match(mean.int_val_list_diffTraits$Var1, df.nlookup$predictor)]
mean.int_val_list_diffTraits$Var2 <- df.nlookup$desc[match(mean.int_val_list_diffTraits$Var2, df.nlookup$predictor)]

#plotting part
# To prepare the long-format data:
if (cohort == "C" | cohort == "D") {
  shap_long <- shap.prep(shap_contrib = data.frame(mean.shap_val[,-ind_bias][ind_subsamp,]),
                         X_train = dataX.rename[ind_subsamp,])
} else {
  shap_long <- shap.prep(shap_contrib = data.frame(mean.shap_val[,-ind_bias]),
                         X_train = dataX.rename)
}

# a plot of dependence with color labelling based on an interacting var
#install_github("kassambara/easyGgplot2")
library(ggplot2)
library(easyGgplot2)

N_graphs <- 6

p <- list()
for (i in 1:N_graphs) {
  p[[i]] <- shap.plot.dependence(data_long = shap_long, 
                                 x = mean.int_val_list_diffTraits$Var1[i], 
                                 color_feature = mean.int_val_list_diffTraits$Var2[i])
}

# Save Dependence plot to file
setwd(paste0(out_path, "/shapley"))
pdf(paste0("Mean_shap_dependence_cohort", cohort, "_plot_top6.pdf"),
    width=7, height=9.5, paper="a4")
ggplot2.multiplot(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]], p[[6]], cols=2)
dev.off()



###
# Plot the Shapley graphs of top 15 and 30 features ranked by p-value
###

# Load the pre-calculated mean absolute shapley values for all features
df.meanAbsShap.obs.out <- read.csv(paste0(out_path, "/pvalue/Mean_shap_cohort", cohort, "_pval.csv"), header = TRUE, stringsAsFactors=FALSE, quote="\"")

# Sort the features by p-values
meanAbsShap.obs.sortby_pval <- df.meanAbsShap.obs.out[,order(df.meanAbsShap.obs.out[2,])]

# Replace the covarite name by full description
colnames(meanAbsShap.obs.sortby_pval) <- df.nlookup$desc[match(colnames(meanAbsShap.obs.sortby_pval), df.nlookup$predictor)]

# SHAP VarImp Plot
pdf(paste0(out_path, "/shapley/Mean_shap_cohort", cohort, "_varimp_plot_top15_by_pval.pdf"), 
    width=7, height=4.75, paper="a4")
xgb.plot.shap(data=dataX.rename,
              shap_contrib=mean.shap_val,
              features=colnames(meanAbsShap.obs.sortby_pval)[1:15],
			                  n_col=3)
pdf(paste0(out_path, "/shapley/Mean_shap_cohort", cohort, "_varimp_plot_top30_by_pval.pdf"),
    width=7, height=9.5, paper="a4")
xgb.plot.shap(data=dataX.rename,
              shap_contrib=mean.shap_val,
              features=colnames(meanAbsShap.obs.sortby_pval)[1:30],
              n_col=3)
dev.off()

