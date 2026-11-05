# Presentation only. Reads saved correlation estimates; never fits a model,
# recomputes a correlation, generates imputations or changes any pupil values.
rebuild_pd_presentation <- function(project_root=getwd(),
    output_root=file.path(project_root,'output'), code_dir=file.path(project_root,'code')) {
  stopifnot(requireNamespace('ggplot2',quietly=TRUE),requireNamespace('ragg',quietly=TRUE))
  project_root <- normalizePath(project_root,mustWork=TRUE)
  spec_env <- new.env(parent=globalenv())
  sys.source(file.path(code_dir,'supplementary_pd_correlations.R'),envir=spec_env)
  spec <- spec_env$pd_correlation_spec()
  stopifnot(unname(spec$pd['RodDiff'])=='EP Net Pupil Change')
  dir.create(file.path(output_root,'figures'),recursive=TRUE,showWarnings=FALSE)
  suppressPackageStartupMessages(library(ggplot2))

  # S2 is a schematic, not an observed participant trace.
  signal <- function(t) -0.2 + 0.2/(1+exp(-(t-1000)/300))
  d <- data.frame(time=seq(-1000,5000,by=2))
  d$pupil <- signal(d$time)
  response <- 4000
  p_onset <- signal(0)
  p_response <- signal(response)
  labels <- c('Early Mean Pupil Size','Early Peak Pupil Change',
              'Late Sustained Dilation','Dilation at Response')
  windows <- data.frame(measure=factor(labels,levels=labels),
    start=c(500,500,2128,response-250),end=c(850,2000,4138,response+250))
  windows <- windows[c(2,1,3,4),]
  fills <- setNames(c('#B5AED6','#C6E3F1','#C7EFCB','#F5D1DA'),labels)
  schematic <- ggplot(d,aes(time,pupil))+
    geom_rect(data=windows,aes(xmin=start,xmax=end,ymin=-Inf,ymax=Inf,fill=measure),
              inherit.aes=FALSE,alpha=0.40,colour=NA)+
    geom_line(aes(colour='Z-scored pupil signal'),linewidth=0.85)+
    geom_vline(aes(xintercept=0,linetype='Stimulus onset'),
               colour='#747474',linewidth=0.7,key_glyph=draw_key_path)+
    geom_vline(aes(xintercept=response,linetype='Participant response'),
               colour='#C43A32',linewidth=0.7,key_glyph=draw_key_path)+
    annotate('segment',x=0,xend=response,y=p_onset,yend=p_onset,
             colour='#B62F28',linetype='dotted',linewidth=0.6)+
    annotate('segment',x=response,xend=response,y=p_onset,yend=p_response,
             colour='#B62F28',linewidth=1.1,
             arrow=grid::arrow(length=grid::unit(0.12,'inches'),type='closed'))+
    annotate('point',x=c(0,response),y=c(p_onset,p_response),colour='#B62F28',size=2.1)+
    annotate('text',x=3700,y=-0.103,label='Net Pupil Change',hjust=1,
             colour='#A92C25',family='Arial',size=4.3)+
    scale_colour_manual(values=c('Z-scored pupil signal'='black'),name=NULL)+
    scale_linetype_manual(values=c('Stimulus onset'='dashed','Participant response'='dashed'),
                          breaks=c('Stimulus onset','Participant response'),name=NULL)+
    scale_fill_manual(values=fills,breaks=labels,name=NULL)+
    scale_x_continuous(breaks=seq(-1000,5000,1000),expand=expansion(mult=0))+
    scale_y_continuous(breaks=c(-0.2,-0.1,0,0.1),expand=expansion(mult=0))+
    coord_cartesian(xlim=c(-1000,5000),ylim=c(-0.245,0.10),clip='off')+
    labs(x='Time from stimulus onset (ms)',y='Z-scored pupil size')+
    guides(colour=guide_legend(order=1),
           linetype=guide_legend(order=2,override.aes=list(colour=c('#747474','#C43A32'))),
           fill=guide_legend(order=3,override.aes=list(alpha=0.40)))+
    theme_classic(base_size=14,base_family='Arial')+
    theme(panel.grid.major=element_line(colour='#E2E2E2',linetype='dotted',linewidth=0.4),
          axis.text=element_text(colour='#303030',size=12.5),axis.title=element_text(size=14),
          legend.position='right',legend.text=element_text(size=12.5),
          legend.key.width=grid::unit(0.40,'inches'),legend.key.height=grid::unit(0.24,'inches'),
          legend.key.spacing.y=grid::unit(0.12,'inches'),legend.spacing.y=grid::unit(0.10,'inches'),
          legend.margin=margin(4,0,4,6),plot.margin=margin(15,15,10,10))
  ggsave(file.path(output_root,'figures/S2_Fig.png'),schematic,
         width=13.6,height=5,units='in',dpi=300,device=ragg::agg_png,bg='white')
  ggsave(file.path(output_root,'figures/S2_Fig.tif'),schematic,
         width=13.6,height=5,units='in',dpi=300,device=ragg::agg_tiff,
         compression='lzw',bg='white')

  # Replot the existing numeric CSV estimates, using the current display labels.
  for(folder in c('correlations_observed','correlations_mi_sensitivity')) {
    results <- utils::read.csv(file.path(project_root,'output',folder,'pd_correlations_all_90.csv'),
                              stringsAsFactors=FALSE,check.names=FALSE)
    stopifnot(nrow(results)==90L)
    out <- file.path(output_root,folder)
    dir.create(out,recursive=TRUE,showWarnings=FALSE)
    for(period_name in c('early_pregnancy','late_pregnancy','postpartum')) {
      x <- results[results$period==period_name,,drop=FALSE]
      stopifnot(nrow(x)==30L)
      p_text <- ifelse(x$p<.001,'<.001',sprintf('%.3f',x$p))
      n_text <- if(x$method[1]=='observed_pairs') paste0('n=',x$n_observed_pairs) else paste0('observed n=',x$n_observed_pairs)
      x$label <- paste0(sprintf('r=%.3f',x$r),'\np=',p_text,'\n',n_text)
      x$font <- ifelse(x$p<=.10,'bold','plain')
      x$pd_label <- factor(unname(spec$pd[x$pd_variable]),levels=rev(unname(spec$pd)))
      symptom_vars <- names(spec$symptoms)[spec$periods==period_name]
      x$symptom_label <- factor(unname(spec$symptoms[x$symptom_variable]),levels=unname(spec$symptoms[symptom_vars]))
      caption <- if(x$method[1]=='observed_pairs')
        'Pearson correlations using observed pairs. Bold: unadjusted p ≤ .10 (exploratory).\nSample size is shown for each pair; 95% intervals and BH-adjusted p-values are in the CSV.' else
        'Correlations pooled across separate full-cohort imputations using Fisher z and Rubin’s rules.\nBold: unadjusted pooled p ≤ .10 (exploratory). Observed pair counts are shown.'
      plot <- ggplot(x,aes(x=symptom_label,y=pd_label,fill=r))+
        geom_tile(color='white')+geom_text(aes(label=label,fontface=font),size=3.1)+
        scale_fill_gradient2(low='#2166AC',mid='white',high='#B2182B',midpoint=0,limits=c(-1,1),name='Pearson r')+
        labs(x=NULL,y=NULL,caption=caption)+theme_minimal(base_size=11)+
        theme(panel.grid=element_blank(),axis.text.x=element_text(face='bold'),
              plot.caption=element_text(hjust=0,size=9))
      ggsave(file.path(out,paste0('pd_correlations_',period_name,'.png')),plot,
             width=10,height=9,dpi=300,bg='white')
      if(folder=='correlations_observed' && period_name=='postpartum') {
        ggsave(file.path(output_root,'figures/S3_Fig.png'),plot,width=10,height=9,dpi=300,bg='white')
        ggsave(file.path(output_root,'figures/S3_Fig.tif'),plot,width=10,height=9,dpi=300,
               device=ragg::agg_tiff,compression='lzw',bg='white')
      }
    }
  }
  invisible(TRUE)
}
