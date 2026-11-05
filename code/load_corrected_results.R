# Safe to source in an existing RStudio session; no imputation or model runs.
.local_source_path <- tryCatch(sys.frame(1)$ofile,error=function(e)NULL)
.local_project_root <- if(!is.null(.local_source_path))normalizePath(file.path(dirname(.local_source_path),'..'))else normalizePath(getwd())
source(file.path(.local_project_root,'code/verify_corrected_project.R'))
verify_corrected_project(.local_project_root)
.local_objects <- readRDS(file.path(.local_project_root,'output/derived/corrected_workspace_objects.rds'))
# Presentation-only rename in the loaded table; saved model objects remain unchanged.
if('master_pd_stability_table' %in% names(.local_objects)) {
  .local_table <- .local_objects$master_pd_stability_table
  .local_table$pd_measure[.local_table$pd_measure=='Rate of Dilation'] <- 'Net Pupil Change'
  .local_objects$master_pd_stability_table <- .local_table
  rm(.local_table)
}
# Replace only this manuscript's old analysis objects and obsolete single-imputed data.
.local_old_names <- unique(c(grep('^(epds_|stais_|stait_|master_)',ls(.GlobalEnv),value=TRUE),
  intersect(c('df_ALL_complete','df_pred_complete','df_gb1_complete','df_gb2_complete',
    'cor_matrix','p_values_matrix','combined_df','group1_data','group2_data',
    'boot_246','boot_stais2','checkpoint_check','imp_246','imp_stais2_baseline','mice_test_epds'),ls(.GlobalEnv))))
if(length(.local_old_names))rm(list=.local_old_names,envir=.GlobalEnv)
invisible(list2env(.local_objects,envir=.GlobalEnv))
assign('raw_file',file.path(.local_project_root,'data/raw/3PADCombined_R_allSubs_25Oct26.xlsx'),envir=.GlobalEnv)
assign('bootstrap_root',file.path(.local_project_root,'output/bootstrap'),envir=.GlobalEnv)
for(.local_key in c('epds','stais','stait'))assign(paste0(.local_key,'_bootstrap_output'),file.path(.local_project_root,'output/bootstrap',.local_key),envir=.GlobalEnv)
rm(.local_key,envir=.GlobalEnv)
source(file.path(.local_project_root,'code/corrected_model_functions.R'),local=.GlobalEnv)
source(file.path(.local_project_root,'code/supplementary_pd_correlations.R'),local=.GlobalEnv)
.local_correlations <- readRDS(file.path(.local_project_root,'output/derived/supplementary_correlations.rds'))
assign('corr_observed',.local_correlations$observed,envir=.GlobalEnv)
assign('corr_mi',.local_correlations$mi,envir=.GlobalEnv)
assign('correlation_sensitivity_comparison',.local_correlations$comparison,envir=.GlobalEnv)
assign('corrected_results_loaded',TRUE,envir=.GlobalEnv)
assign('corrected_results_project',.local_project_root,envir=.GlobalEnv)
message('Verified corrected results loaded: 1,000 paired bootstraps per outcome, corrected tables/figures, observed correlations and MI sensitivity.')
rm(list=c('.local_source_path','.local_project_root','.local_objects','.local_old_names','.local_correlations'),envir=.GlobalEnv)
