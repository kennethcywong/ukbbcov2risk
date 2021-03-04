#!/usr/bin/env Rscript

library("argparse")
suppressPackageStartupMessages(library("argparse"))

help_message <-
  "usage: Rscript analysis_dataset.R
        <--cohort A>
        <--dspath /home/kenneth/out/UKBB/covid19/run_03Nov2020/impute_with_replace/xgb_model>
        <--outpath /home/kenneth/out/UKBB/covid19/run_03Nov2020/impute_with_replace>
\nRequired:
        --cohort    Cohort (A,B,C,D) to be analysed
        --dspath    Folder storing training data and XGboost built model
        --outpath   Output folder
        --help | -h Display this help message\n"

parser <- ArgumentParser()
parser$add_argument("--cohort", action="store", type="character", default="A", dest="cohort")
parser$add_argument("--dspath", action="store")
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

#assumed to be followed by fold number and .Rdata
N_folds = 5 
cohort = args$cohort

dataX_path = paste0(args$dspath, "/dataX_trainTestcomb_cohort", cohort, ".Rdata")
xgbmodel_path = paste0(args$dspath, "/xgb_model_cohort_", cohort, "_fold_")
out_path = args$outpath # "/home/kenneth/out/UKBB/covid19/run_03Nov2020/impute_with_replace"

# Load the predictor-description lookup table
df.nlookup <- read.csv("/home/kenneth/data/UKBB/covid19/meta/names_lookup.csv", header = TRUE, stringsAsFactors=FALSE, quote="\"")

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
mean_list.varimp$Feature <- df.nlookup$desc[match(mean_list.varimp$Feature, df.nlookup$predictor)]
mean_list.varimp  = arrange(mean_list.varimp, -Gain)
library(data.table)
mean_list.varimp = as.data.table(mean_list.varimp)


#_______________________________________________
#**Mean XGboost VarImp plot**
#______________________________________________
library(Ckmeans.1d.dp)
library(ggplot2)

# Set the output path of saved graph and Rdata
setwd(paste0(out_path, "/xgb_model"))

# Plot the varImp (top 30 features)
#pdf(paste0("Mean_xgb.importance_cohort", cohort, "_varimp_plot_top30_gg", ".pdf"),
#    width=7, height=7.5)
p <- xgb.ggplot.importance(importance_matrix=mean_list.varimp[1:30])
p+theme_bw()+
  theme(axis.text=element_text(face="plain", angle=0, size="11", color="black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none")+
  labs(x="Risk Factor", y="XGBoost VarImp", title="")
ggsave(device="png", file=paste0("Mean_xgb.importance_cohort", cohort, "_varimp_plot_top30", ".png"), width=6, height=7, unit="in", dpi=300)
dev.off()

# Plot the varImp (top 15 features)
p <- xgb.ggplot.importance(importance_matrix=mean_list.varimp[1:15])
p+theme_bw()+
  theme(axis.text=element_text(face="plain", angle=0, size="10", color="black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none")+
  labs(x="Risk Factor", y="XGBoost VarImp", title="")
ggsave(device="png", file=paste0("Mean_xgb.importance_cohort", cohort, "_varimp_plot_top15", ".png"), width=4, height=3.5, unit="in", dpi=300)
dev.off()

# Plot the varImp (top 15 features)
p <- xgb.ggplot.importance(importance_matrix=mean_list.varimp[1:10])
p+theme_bw()+
  theme(axis.text=element_text(face="plain", angle=0, size="10", color="black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none")+
  labs(x="Risk Factor", y="XGBoost VarImp", title="")
ggsave(device="png", file=paste0("Mean_xgb.importance_cohort", cohort, "_varimp_plot_top10", ".png"), width=4, height=2.2, unit="in", dpi=300)
dev.off()

# Plot the varImp (top 5 features)
p <- xgb.ggplot.importance(importance_matrix=mean_list.varimp[1:5])
p+theme_bw()+
  theme(axis.text=element_text(face="plain", angle=0, size="10", color="black"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none")+
  labs(x="Risk Factor", y="XGBoost VarImp", title="")
ggsave(device="png", file=paste0("Mean_xgb.importance_cohort", cohort, "_varimp_plot_top5", ".png"), width=4, height=1.7, unit="in", dpi=300)
dev.off()

#png(paste0("Mean_xgb.importance_cohort", cohort, "_varimp_plot_top15", ".png"), units="in", res=600
#    ,width=8, height=8) #, paper="a4")
#xgb.plot.importance(importance_matrix=mean_list.varimp[1:15], cex.lab=1.5, cex.axis=1.5)
#	xlab="XGBoost VarImp", ylab="Clinical Risk", col="blue", cex.lab=1.5)
#dev.off()



