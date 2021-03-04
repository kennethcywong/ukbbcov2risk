#!/bin/bash

VER="14Dec2020"
RUN="30Dec2020"

Rscript analysis_dataset.R --pval TRUE --cohort A --dspath /home/kenneth/out/UKBB/covid19/rerun_${RUN}/df_cohorts_${VER}.RData --metapath /home/kenneth/data/UKBB/covid19/meta/${RUN} --cpu 32 --outpath /home/kenneth/out/UKBB/covid19/rerun_${RUN} > /home/kenneth/out/UKBB/covid19/rerun_${RUN}/cohortA.log 2>&1
Rscript analysis_dataset.R --pval TRUE --cohort B --dspath /home/kenneth/out/UKBB/covid19/rerun_${RUN}/df_cohorts_${VER}.RData --metapath /home/kenneth/data/UKBB/covid19/meta/${RUN} --cpu 32 --outpath /home/kenneth/out/UKBB/covid19/rerun_${RUN} > /home/kenneth/out/UKBB/covid19/rerun_${RUN}/cohortB.log 2>&1
Rscript analysis_dataset.R --pval TRUE --cohort C --dspath /home/kenneth/out/UKBB/covid19/rerun_${RUN}/df_cohorts_${VER}.RData --metapath /home/kenneth/data/UKBB/covid19/meta/${RUN} --cpu 32 --outpath /home/kenneth/out/UKBB/covid19/rerun_${RUN} > /home/kenneth/out/UKBB/covid19/rerun_${RUN}/cohortC.log 2>&1
Rscript analysis_dataset.R --pval TRUE --cohort D --dspath /home/kenneth/out/UKBB/covid19/rerun_${RUN}/df_cohorts_${VER}.RData --metapath /home/kenneth/data/UKBB/covid19/meta/${RUN} --cpu 32--outpath /home/kenneth/out/UKBB/covid19/rerun_${RUN} > /home/kenneth/out/UKBB/covid19/rerun_${RUN}/cohortD.log 2>&1

Rscript analysis_shapley.R --cohort A --dspath /home/kenneth/out/UKBB/covid19/rerun_${RUN}/xgb_model --outpath /home/kenneth/out/UKBB/covid19/rerun_${RUN} >> /home/kenneth/out/UKBB/covid19/rerun_${RUN}/cohortA.log 2>&1
Rscript analysis_shapley.R --cohort B --dspath /home/kenneth/out/UKBB/covid19/rerun_${RUN}/xgb_model --outpath /home/kenneth/out/UKBB/covid19/rerun_${RUN} >> /home/kenneth/out/UKBB/covid19/rerun_${RUN}/cohortB.log 2>&1
Rscript analysis_shapley.R --cohort C --dspath /home/kenneth/out/UKBB/covid19/rerun_${RUN}/xgb_model --outpath /home/kenneth/out/UKBB/covid19/rerun_${RUN} >> /home/kenneth/out/UKBB/covid19/rerun_${RUN}/cohortC.log 2>&1
Rscript analysis_shapley.R --cohort D --dspath /home/kenneth/out/UKBB/covid19/rerun_${RUN}/xgb_model --outpath /home/kenneth/out/UKBB/covid19/rerun_${RUN} >> /home/kenneth/out/UKBB/covid19/rerun_${RUN}/cohortD.log 2>&1
