additional_library <- Sys.getenv("THREEPAD_ADDITIONAL_R_LIB")
if (nzchar(additional_library)) .libPaths(c(.libPaths(), additional_library))
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) >= 3L)
root <- normalizePath(args[1])
project <- normalizePath(args[2])
corrected_file <- normalizePath(args[3])
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])
refit_code <- dirname(normalizePath(script))
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(purrr);library(tibble);library(mice);library(glmnet);library(pROC)})
source(file.path(refit_code,'original_workflow_functions.R'));data<-readRDS(corrected_file)
key<-args[4];b<-if(key=='epds')220L else 295L;outcome<-if(key=='epds')'outcome_epds' else 'outcome_stai_s'

base_fit<-fit_elastic_net
repairs<-list();j<-0L
fit_elastic_net<-function(...) {
 j<<-j+1L;args<-list(...);old<-base_fit(...)
 if(!grepl('Convergence|convergence issue',old$glmnet_warning_text))return(old)
 matrices<-make_model_matrices(args$train_data,args$validation_data,args$predictor_vars)
 y<-args$train_data[[args$outcome_var]]
 warnings<-character()
 capture<-function(expr)withCallingHandlers(expr,warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
 path<-capture(glmnet::glmnet(matrices$train,y,family='binomial',alpha=old$alpha,standardize=TRUE))
 lambdas<-sort(unique(c(path$lambda[path$lambda>=old$lambda],old$lambda)),decreasing=TRUE)
 fit<-capture(glmnet::glmnet(matrices$train,y,family='binomial',alpha=old$alpha,lambda=lambdas,standardize=TRUE))
 stopifnot(fit$jerr==0L,any(abs(fit$lambda-old$lambda)<1e-14))
 new<-old
 new$train_prediction<-as.numeric(predict(fit,newx=matrices$train,s=old$lambda,type='response'))
 new$validation_prediction<-as.numeric(predict(fit,newx=matrices$validation,s=old$lambda,type='response'))
 co<-as.matrix(coef(fit,s=old$lambda));new$coefficients<-tibble(variable=rownames(co),coefficient=co[,1]) %>% filter(variable!='(Intercept)') %>% mutate(selected=coefficient!=0,direction=case_when(coefficient>0~'positive',coefficient<0~'negative',TRUE~'zero'))
 # Cross-validation diagnostics are retained; repair details are recorded separately.
 repairs[[length(repairs)+1L]]<<-list(imputation=j,alpha=old$alpha,lambda=old$lambda,original_selected=sum(old$coefficients$selected),repaired_selected=sum(new$coefficients$selected),glmnet_jerr=fit$jerr,repair_warnings=warnings,original=old,repaired=new)
 cat(key,b,'imputation',j,'path refit jerr',fit$jerr,'selected',sum(new$coefficients$selected),'\n');flush.console();new
}
full<-run_one_bootstrap(data,outcome,full_vars,m=10,bootstrap_id=b,seed=26092026)
saved<-readRDS(file.path(root,'bootstrap',key,sprintf('%s_bootstrap_%04d.rds',key,b)))
stopifnot(identical(full$bootstrap_index,saved$bootstrap_index))
base<-lapply(saved[c('performance','tuning','coefficients','diagnostics')],function(x)x%>%filter(model=='baseline')%>%select(-model));base$predictions<-saved$predictions$baseline;base$bootstrap_index<-saved$bootstrap_index
run_one_bootstrap<-function(data,outcome_var,predictor_vars,...)if(identical(predictor_vars,baseline_vars))base else full
paired<-run_one_paired_bootstrap(data,outcome,baseline_vars,full_vars,m=10,bootstrap_id=b,seed=26092026)
repairdir<-file.path(root,'numerical_repairs',key);dir.create(repairdir,recursive=TRUE,showWarnings=FALSE)
saveRDS(paired,file.path(repairdir,sprintf('%s_bootstrap_%04d.rds',key,b)))
saveRDS(repairs,file.path(repairdir,'repair_details.rds'))
cat('REPAIR COMPLETE',key,'\n')
