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
stopifnot(requireNamespace('brglm2',quietly=TRUE))

source(file.path(refit_code,'original_workflow_functions.R'))
original <- readRDS(file.path(project,'data/derived/analysis_unimputed.rds'))
corrected <- readRDS(corrected_file)
stopifnot(identical(original[setdiff(names(original),c('PARDiff','PARDiff2','RodDiff','RodDiff2'))],corrected[setdiff(names(corrected),c('PARDiff','PARDiff2','RodDiff','RodDiff2'))]))
original_run_one_bootstrap <- run_one_bootstrap
extract_baseline <- function(saved) {
 out <- lapply(saved[c('performance','tuning','coefficients','diagnostics')],function(x) x %>% filter(model=='baseline') %>% select(-model))
 out$predictions <- saved$predictions$baseline; out$bootstrap_index <- saved$bootstrap_index
 out
}
run_one_bootstrap <- function(data,outcome_var,predictor_vars,m=10,bootstrap_id=1,seed=26092026) {
 if (identical(predictor_vars,baseline_vars)) {
  stopifnot(identical(data[baseline_vars],original[baseline_vars]))
  return(extract_baseline(current_saved))
 }
 original_run_one_bootstrap(data,outcome_var,predictor_vars,m,bootstrap_id,seed)
}
args <- args[-seq_len(3L)]
mode <- if(length(args)) args[1] else 'check'
keys <- c(epds='outcome_epds',stais='outcome_stai_s',stait='outcome_stai_t')
if (mode=='check') {
 writeLines(capture.output(sessionInfo()),file.path(root,'sessionInfo.txt'))
 saved <- readRDS(file.path(project,'output/bootstrap/epds/epds_bootstrap_0001.rds'))
 t <- proc.time()[3]
 reproduced <- original_run_one_bootstrap(original,keys[['epds']],full_vars,m=10,bootstrap_id=1,seed=26092026)
 stopifnot(isTRUE(all.equal(reproduced$predictions,saved$predictions$full,tolerance=1e-12)))
 stopifnot(isTRUE(all.equal(reproduced$performance,saved$performance %>% filter(model=='full') %>% select(-model),tolerance=1e-12)))
 saveRDS(list(passed=TRUE,elapsed_seconds=unname(proc.time()[3]-t),glmnet=as.character(packageVersion('glmnet')),mice=as.character(packageVersion('mice')),brglm2=as.character(packageVersion('brglm2'))),file.path(root,'original_reproduction_check.rds'))
 cat('PASS: original full EPDS bootstrap reproduced, including predictions and calibration.\n')
 quit(save='no')
}
stopifnot(file.exists(file.path(root,'original_reproduction_check.rds')))
key <- args[2]; stopifnot(key %in% names(keys)); outcome <- keys[[key]]
if(mode=='pilot') ids <- 1L else if(mode=='run') ids <- seq_len(1000L) else stop('Unknown mode')
outdir <- file.path(root,'bootstrap',key); dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
for(b in ids) {
 target <- file.path(outdir,sprintf('%s_bootstrap_%04d.rds',key,b))
 if(file.exists(target)) { tmp <- readRDS(target); stopifnot(tmp$comparison$bootstrap==b); next }
 current_saved <- readRDS(file.path(project,'output/bootstrap',key,sprintf('%s_bootstrap_%04d.rds',key,b)))
 t <- proc.time()[3]
 result <- run_one_paired_bootstrap(corrected,outcome,baseline_vars,full_vars,m=10,bootstrap_id=b,seed=26092026)
 stopifnot(identical(result$predictions$baseline,current_saved$predictions$baseline))
 stopifnot(identical(result$bootstrap_index,current_saved$bootstrap_index))
 stopifnot(!anyNA(result$performance$auc_bootstrap),!anyNA(result$performance$auc_validation))
 temporary <- paste0(target,'.tmp'); saveRDS(result,temporary); stopifnot(file.rename(temporary,target))
 cat(format(Sys.time(),'%Y-%m-%d %H:%M:%S'),key,b,'saved;',round(proc.time()[3]-t,1),'seconds\n'); flush.console()
}
apparent_path <- file.path(root,paste0(key,'_apparent_full_corrected.rds'))
if(mode=='run' && !file.exists(apparent_path)) {
 apparent <- fit_apparent_model(corrected,outcome,full_vars,m=10,seed=26092026)
 saveRDS(apparent,apparent_path)
}
cat('COMPLETE',mode,key,'\n')
