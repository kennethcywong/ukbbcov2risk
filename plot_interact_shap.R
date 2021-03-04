cohort <- "C"
N_folds <- 5

setwd("/home/kenneth/data/UKBB/covid19/")
dspath <- "/home/kenneth/out/UKBB/covid19/rerun_30Dec2020"
outpath <- "/home/kenneth/out/UKBB/covid19/rerun_tmp"
xgbmodel_path <- paste0(outpath, "/xgb_model/xgb_model_cohort_", cohort, "_fold_")

df.nlookup <- read.csv("meta/names_lookup_dot.csv", header = TRUE, stringsAsFactors=FALSE, quote="\"")

load(paste0(outpath, "/xgb_model/dataX_trainTestcomb_cohort", cohort, ".Rdata"))

no_features = ncol(dataX)
feature_names = colnames(dataX)
xgb.models = list()

for (i in 1:N_folds) {
  load(paste0(xgbmodel_path, i, ".Rdata"))
  xgb.models[[i]]  <-  best_xgb_model
}

library(SHAPforxgboost)
library(xgboost)
library(dplyr)
setwd(paste0(outpath, "/shapley"))

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

#** sort the variables by mean |shap| **
ind_bias = which(colnames(mean.shap_val) == "BIAS")
mat_meanAbsShap = data.frame( Feature = colnames(mean.shap_val[,-ind_bias]) ,
                              meanAbsShap = colMeans( abs(mean.shap_val[,-ind_bias]) )
			   )


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
mean.int_val_list_diffTraits =    arrange(mean.int_val_list[ind_interact,], -value)

setwd(paste0(outpath, "/shapley"))
save(ind_subsamp, mean.int_val,mean.int_val_list_diffTraits, file= paste0("Interaction_Shap", cohort, ".Rdata"))


setwd(paste0(dspath, "/shapley"))

if (cohort == "A" | cohort == "B") {
  load(paste0("Interaction_Shap", cohort, ".Rdata"))
  mean.int_val_list_diffTraits$Var1 = as.character(mean.int_val_list_diffTraits$Var1)
  mean.int_val_list_diffTraits$Var2 = as.character(mean.int_val_list_diffTraits$Var2)
}

# Lookup the description name
mean.int_val_list_diffTraits$Var1 <- df.nlookup$desc[match(mean.int_val_list_diffTraits$Var1, df.nlookup$predictor)]
mean.int_val_list_diffTraits$Var2 <- df.nlookup$desc[match(mean.int_val_list_diffTraits$Var2, df.nlookup$predictor)]


colnames(mean.shap_val)[-ind_bias] <- df.nlookup$desc[match(colnames(mean.shap_val[,-ind_bias]), df.nlookup$predictor)]

colnames(dataX) <- df.nlookup$desc[match(colnames(dataX), df.nlookup$predictor)]

# To prepare the long-format data:
if (cohort == "C" | cohort == "D") {
  shap_long <- shap.prep(shap_contrib = data.frame(mean.shap_val[,-ind_bias][ind_subsamp,]),
				                          X_train = dataX[ind_subsamp,])
} else {
  shap_long <- shap.prep(shap_contrib = data.frame(mean.shap_val[,-ind_bias]),
				                          X_train = dataX)
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
setwd(paste0(outpath, "/shapley"))
pdf(paste0("Mean_shap_dependence_cohort", cohort, "_plot_top6.pdf"),
        width=7, height=9.5, paper="a4")
ggplot2.multiplot(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]],
		                    p[[6]], cols=2)
dev.off()

