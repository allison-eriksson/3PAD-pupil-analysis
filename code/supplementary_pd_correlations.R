# Supplementary PD correlations, independent of elastic-net/bootstrap fitting.
# Source this file once. No analysis runs until run_pd_correlations() is called.
# Methods must be selected explicitly: observed_pairs or multiple_imputation.

pd_correlation_spec <- function() {
  pd <- c(PTRDiff='EP Early Mean Pupil Size', PTRDiff2='LP Early Mean Pupil Size',
    PCUDiff='EP Early Peak Pupil Change', PCUDiff2='LP Early Peak Pupil Change',
    PAADiff='EP Late Sustained Dilation', PAADiff2='LP Late Sustained Dilation',
    PARDiff='EP Dilation at Response', PARDiff2='LP Dilation at Response',
    RodDiff='EP Net Pupil Change', RodDiff2='LP Net Pupil Change')
  symptoms <- c(gb1_epds_sum='EP EPDS',gb1_stai_s_sum='EP STAI-S',gb1_stai_t_sum='EP STAI-T',
    gb2_epds_sum='LP EPDS',gb2_stai_s_sum='LP STAI-S',gb2_stai_t_sum='LP STAI-T',
    ppv6_epds_sum='Postpartum EPDS',ppv6_stai_s_sum='Postpartum STAI-S',ppv6_stai_t_sum='Postpartum STAI-T')
  periods <- rep(c('early_pregnancy','late_pregnancy','postpartum'),each=3)
  names(periods) <- names(symptoms)
  list(pd=pd,symptoms=symptoms,periods=periods)
}

prepare_pd_correlation_data <- function(raw_data, auxiliary_vars=NULL) {
  spec <- pd_correlation_spec()
  needed <- unique(c('pid',names(spec$pd),names(spec$symptoms),auxiliary_vars))
  stopifnot(all(needed %in% names(raw_data)),!anyDuplicated(raw_data$pid))
  # The same cohort as the prediction analyses; retain row-wise alignment.
  cohort <- raw_data[stats::complete.cases(raw_data[c('ppv6_epds_sum',
    'ppv6_stai_s_sum','ppv6_stai_t_sum')]),needed,drop=FALSE]
  stopifnot(nrow(cohort)==91L)
  for(v in c(names(spec$pd),names(spec$symptoms))) {
    original <- cohort[[v]]
    cohort[[v]] <- as.numeric(original)
    stopifnot(!any(!is.na(original) & is.na(cohort[[v]])))
  }
  for(v in auxiliary_vars) if(is.character(cohort[[v]])) cohort[[v]] <- factor(cohort[[v]])
  cohort
}

observed_pd_correlations <- function(correlation_data) {
  spec <- pd_correlation_spec();rows <- list()
  for(pd_var in names(spec$pd)) for(symptom_var in names(spec$symptoms)) {
    keep <- stats::complete.cases(correlation_data[c(pd_var,symptom_var)])
    n <- sum(keep);x <- correlation_data[[pd_var]][keep];y <- correlation_data[[symptom_var]][keep]
    stopifnot(n>=4L,stats::sd(x)>0,stats::sd(y)>0)
    fit <- stats::cor.test(x,y,method='pearson',conf.level=.95)
    rows[[length(rows)+1L]] <- data.frame(pd_variable=pd_var,symptom_variable=symptom_var,
      period=unname(spec$periods[symptom_var]),method='observed_pairs',n_analysis=n,
      n_observed_pairs=n,r=unname(fit$estimate),lower_95=fit$conf.int[1],
      upper_95=fit$conf.int[2],p=fit$p.value,fmi=NA_real_,stringsAsFactors=FALSE)
  }
  do.call(rbind,rows)
}

make_correlation_imputations <- function(correlation_data,m=50L,maxit=20L,
                                         seed=26102026L,predictor_matrix=NULL) {
  # Separate full-cohort imputations, never the resampled/validation imputations.
  # This transparent starting specification includes the 10 PD indices and
  # all 9 CONTINUOUS symptom scores examined in the 90 correlations.
  # Do not replace these continuous scores with binary clinical thresholds.
  # Optionally supply a reviewed predictor matrix. Review trace/density plots,
  # logged events and plausibility before using MI results for publication.
  stopifnot(requireNamespace('mice',quietly=TRUE),m>=2L,maxit>=1L)
  d <- correlation_data[setdiff(names(correlation_data),'pid')]
  method <- vapply(d,function(x)if(!anyNA(x))'' else if(is.numeric(x))'pmm' else 'cart',character(1))
  if(is.null(predictor_matrix)) {
    predictor_matrix <- mice::make.predictorMatrix(d)
    diag(predictor_matrix) <- 0L
  }
  stopifnot(identical(dim(predictor_matrix),c(ncol(d),ncol(d))))
  imp <- mice::mice(d,m=m,maxit=maxit,method=method,
    predictorMatrix=predictor_matrix,seed=seed,printFlag=FALSE)
  # The postpartum continuous outcomes are fully observed and remain unchanged.
  stopifnot(all(imp$method[c('ppv6_epds_sum','ppv6_stai_s_sum','ppv6_stai_t_sum')]==''))
  imp
}

pooled_pd_correlations <- function(correlation_data,imputations) {
  # Fisher-z estimates and 1/(n-3) complete-data variances, pooled by Rubin's
  # rules using mice::pool.scalar. Rubin t degrees of freedom are used with
  # the default large-sample complete-data reference. This is an approximate
  # Fisher-z inference, not a pooled set of ordinary cor.test() p-values.
  # References: https://stefvanbuuren.name/fimd/sec-pooling.html
  # https://amices.org/mice/reference/pool.scalar.html
  stopifnot(inherits(imputations,'mids'),imputations$m>=2L)
  spec <- pd_correlation_spec();completed <- mice::complete(imputations,'all');rows <- list()
  stopifnot(all(vapply(completed,nrow,integer(1))==nrow(correlation_data)))
  for(pd_var in names(spec$pd)) for(symptom_var in names(spec$symptoms)) {
    r <- vapply(completed,function(d){
      stopifnot(!anyNA(d[c(pd_var,symptom_var)]))
      stats::cor(d[[pd_var]],d[[symptom_var]],method='pearson')
    },numeric(1))
    stopifnot(all(is.finite(r)),all(abs(r)<1))
    n <- nrow(correlation_data)
    pool <- mice::pool.scalar(Q=atanh(r),U=rep(1/(n-3),length(r)))
    se <- sqrt(pool$t);critical <- stats::qt(.975,df=pool$df)
    p <- 2*stats::pt(abs(pool$qbar/se),df=pool$df,lower.tail=FALSE)
    rows[[length(rows)+1L]] <- data.frame(pd_variable=pd_var,symptom_variable=symptom_var,
      period=unname(spec$periods[symptom_var]),method='multiple_imputation',n_analysis=n,
      n_observed_pairs=sum(stats::complete.cases(correlation_data[c(pd_var,symptom_var)])),
      r=tanh(pool$qbar),lower_95=tanh(pool$qbar-critical*se),
      upper_95=tanh(pool$qbar+critical*se),p=p,fmi=pool$fmi,stringsAsFactors=FALSE)
  }
  do.call(rbind,rows)
}

write_pd_correlation_outputs <- function(results,output_dir) {
  stopifnot(requireNamespace('ggplot2',quietly=TRUE))
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  spec <- pd_correlation_spec()
  # Retain raw p-values used by the old supplementary figure. Also export
  # BH-adjusted p-values for each predeclared set of 30 comparisons.
  results$p_BH_within_period <- ave(results$p,results$period,FUN=function(x)stats::p.adjust(x,'BH'))
  utils::write.csv(results,file.path(output_dir,'pd_correlations_all_90.csv'),row.names=FALSE)
  for(period_name in unique(results$period)) {
    d <- results[results$period==period_name,,drop=FALSE]
    utils::write.csv(d,file.path(output_dir,paste0('pd_correlations_',period_name,'.csv')),row.names=FALSE)
    p_text <- ifelse(d$p<.001,'<.001',sprintf('%.3f',d$p))
    n_text <- if(d$method[1]=='observed_pairs')paste0('n=',d$n_observed_pairs) else paste0('observed n=',d$n_observed_pairs)
    d$label <- paste0(sprintf('r=%.3f',d$r),'\np=',p_text,'\n',n_text)
    d$font <- ifelse(d$p<=.10,'bold','plain')
    d$pd_label <- factor(unname(spec$pd[d$pd_variable]),levels=rev(unname(spec$pd)))
    labels <- names(spec$symptoms)[spec$periods==period_name]
    d$symptom_label <- factor(unname(spec$symptoms[d$symptom_variable]),levels=unname(spec$symptoms[labels]))
    caption <- if(d$method[1]=='observed_pairs')
      'Pearson correlations using observed pairs. Bold: unadjusted p ≤ .10 (exploratory).\nSample size is shown for each pair; 95% intervals and BH-adjusted p-values are in the CSV.' else
      'Correlations pooled across separate full-cohort imputations using Fisher z and Rubin’s rules.\nBold: unadjusted pooled p ≤ .10 (exploratory). Observed pair counts are shown.'
    plot <- ggplot2::ggplot(d,ggplot2::aes(x=symptom_label,y=pd_label,fill=r))+
      ggplot2::geom_tile(color='white')+
      ggplot2::geom_text(ggplot2::aes(label=label,fontface=font),size=3.1)+
      ggplot2::scale_fill_gradient2(low='#2166AC',mid='white',high='#B2182B',midpoint=0,limits=c(-1,1),name='Pearson r')+
      ggplot2::labs(x=NULL,y=NULL,caption=caption)+ggplot2::theme_minimal(base_size=11)+
      ggplot2::theme(panel.grid=ggplot2::element_blank(),axis.text.x=ggplot2::element_text(face='bold'),
        plot.caption=ggplot2::element_text(hjust=0,size=9))
    ggplot2::ggsave(file.path(output_dir,paste0('pd_correlations_',period_name,'.png')),plot,
      width=10,height=9,dpi=300,bg='white')
  }
  invisible(results)
}

run_pd_correlations <- function(raw_data,output_dir,method,m=50L,maxit=20L,
                               seed=26102026L,predictor_matrix=NULL,auxiliary_vars=NULL) {
  # Requiring an explicit method avoids changing the missing-data strategy silently.
  stopifnot(length(method)==1L,method %in% c('observed_pairs','multiple_imputation'))
  d <- prepare_pd_correlation_data(raw_data,auxiliary_vars)
  if(method=='observed_pairs') results <- observed_pd_correlations(d) else {
    imp <- make_correlation_imputations(d,m,maxit,seed,predictor_matrix)
    dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
    saveRDS(imp,file.path(output_dir,'correlation_imputations.rds'))
    utils::write.csv(as.data.frame(imp$loggedEvents),file.path(output_dir,'mice_logged_events.csv'),row.names=FALSE)
    results <- pooled_pd_correlations(d,imp)
    if(!is.null(imp$loggedEvents))warning('Inspect MICE logged events before interpreting these MI correlations.')
  }
  write_pd_correlation_outputs(results,output_dir)
}
