additional_library <- Sys.getenv("THREEPAD_ADDITIONAL_R_LIB")
if (nzchar(additional_library)) .libPaths(c(.libPaths(), additional_library))
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) >= 3L)
root <- normalizePath(args[1])
project <- normalizePath(args[2])
corrected_file <- normalizePath(args[3])
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])
refit_code <- dirname(normalizePath(script))
suppressPackageStartupMessages({library(dplyr);library(tibble)})
source(file.path(refit_code,'original_workflow_functions.R'));rows<-list()
for(key in c('epds','stais','stait')) {
 b<-switch(key,epds=220L,stais=295L,stait=616L);name<-sprintf('%s_bootstrap_%04d.rds',key,b)
 old<-readRDS(file.path(root,'bootstrap',key,name));new<-readRDS(file.path(root,'numerical_repairs',key,name))
 stopifnot(identical(old$bootstrap_index,new$bootstrap_index),identical(old$predictions$baseline,new$predictions$baseline),identical(old$performance[old$performance$model=='baseline',],new$performance[new$performance$model=='baseline',]))
 if(key!='stait') {
  repairs<-readRDS(file.path(root,'numerical_repairs',key,'repair_details.rds'));stopifnot(length(repairs)==1,repairs[[1]]$glmnet_jerr==0)
  for(s in c('bootstrap','validation')) {
   pred<-new$predictions$full[[s]]$predicted_prob
   oldpred<-old$predictions$full[[s]]$predicted_prob
   resultname<-if(s=='bootstrap')'train_prediction' else 'validation_prediction'
   expected<-oldpred+(repairs[[1]]$repaired[[resultname]]-repairs[[1]]$original[[resultname]])/10
   stopifnot(isTRUE(all.equal(pred,expected,tolerance=1e-12)))
   idi<-calculate_idi(new$predictions$full[[s]]$outcome,new$predictions$baseline[[s]]$predicted_prob,pred)
   stopifnot(isTRUE(all.equal(idi$idi,new$comparison[[paste0('idi_',s)]],tolerance=1e-12)))
  }
 }
 for(n in c('auc_bootstrap','auc_validation','calibration_slope_bootstrap','calibration_slope_validation')) {
  a<-old$performance[[n]][old$performance$model=='full'];z<-new$performance[[n]][new$performance$model=='full']
  rows[[length(rows)+1L]]<-tibble(outcome=key,bootstrap=b,quantity=n,original=a,repaired=z,difference=z-a)
 }
}
write.csv(bind_rows(rows),file.path(root,'audit/numerical_repair_summary.csv'),row.names=FALSE);cat('All repair pairing, pooled-prediction and IDI checks passed.\n')
