# Run this file alone to reproduce the supplementary analyses; never reruns elastic net.
if(!requireNamespace('here',quietly=TRUE))stop('Install the project dependency here first.')
source(here::here('code','verify_corrected_project.R'))
verify_corrected_project(here::here())
source(here::here('code','supplementary_pd_correlations.R'))
correlation_raw <- readxl::read_excel(here::here('data','raw','3PADCombined_R_allSubs_25Oct26.xlsx'))
corr_observed <- run_pd_correlations(correlation_raw,here::here('output','correlations_observed'),method='observed_pairs')
corr_mi <- run_pd_correlations(correlation_raw,here::here('output','correlations_mi_sensitivity'),
  method='multiple_imputation',m=50L,maxit=20L,seed=26102026L)
correlation_sensitivity_comparison <- merge(corr_observed,corr_mi,
  by=c('pd_variable','symptom_variable','period'),suffixes=c('_observed','_mi'))
correlation_sensitivity_comparison$delta_r <- correlation_sensitivity_comparison$r_mi-correlation_sensitivity_comparison$r_observed
correlation_sensitivity_comparison$nominal_p05_same <- (correlation_sensitivity_comparison$p_observed<.05)==(correlation_sensitivity_comparison$p_mi<.05)
utils::write.csv(correlation_sensitivity_comparison,here::here('output','correlations_comparison.csv'),row.names=FALSE)
saveRDS(list(observed=corr_observed,mi=corr_mi,comparison=correlation_sensitivity_comparison),
  here::here('output','derived','supplementary_correlations.rds'))
source(here::here('code','check_correlation_imputations.R'))
check_correlation_imputations(here::here())
# Refresh hashes only for these explicitly regenerated correlation outputs.
source(here::here('code','refresh_correlation_manifest.R'))
message('Supplementary observed-pair correlations and separate MI sensitivity analysis saved.')
