check_correlation_imputations <- function(project_root) {
  imp<-readRDS(file.path(project_root,'output/correlations_mi_sensitivity/correlation_imputations.rds'))
  cat('Logged events:',if(is.null(imp$loggedEvents))0 else nrow(imp$loggedEvents),'\n')
  conv<-mice::convergence(imp,diagnostic="ac")
  write.csv(conv,file.path(project_root,'output/correlations_mi_sensitivity/convergence_diagnostics.csv'),row.names=FALSE)
  diagnostic_summary<-list()
  for(v in names(imp$method)[imp$method!='']) {
   png(file.path(project_root,'output/correlations_mi_sensitivity',paste0('trace_',v,'.png')),width=1200,height=700)
   print(plot(imp,c(v)));dev.off()
   observed_values<-imp$data[[v]][!is.na(imp$data[[v]])]
   imputed_values<-as.numeric(unlist(imp$imp[[v]]))
   stopifnot(all(is.finite(imputed_values)),all(imputed_values %in% observed_values))
   density_data<-data.frame(value=c(observed_values,imputed_values),source=c(rep('Observed',length(observed_values)),rep('Imputed (pooled)',length(imputed_values))))
   p<-ggplot(density_data,aes(value,color=source))+geom_density(linewidth=.8)+theme_classic()+labs(x=v,y='Density',color=NULL)
   ggsave(file.path(project_root,'output/correlations_mi_sensitivity',paste0('density_',v,'.png')),p,width=8,height=5,dpi=130)
   means<-imp$chainMean[v,,]
   diagnostic_summary[[v]]<-data.frame(variable=v,n_missing=nrow(imp$imp[[v]]),observed_min=min(observed_values),observed_max=max(observed_values),imputed_min=min(imputed_values),imputed_max=max(imputed_values),all_imputed_values_are_observed_donors=TRUE,mean_first_10=mean(means[1:10,]),mean_last_10=mean(means[11:20,]),mean_chain_sd=mean(apply(means,2,sd)))
  }
  write.csv(do.call(rbind,diagnostic_summary),file.path(project_root,'output/correlations_mi_sensitivity/imputation_diagnostics_summary.csv'),row.names=FALSE)
  
  invisible(TRUE)
}
