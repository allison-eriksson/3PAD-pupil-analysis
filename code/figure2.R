library(ggplot2)
d<-read.csv(file.path(p,'master_calibration_plot_data.csv'))
l<-read.csv(file.path(p,'master_calibration_labels.csv'))
for(x in c('d','l')){v<-get(x);v$outcome<-factor(v$outcome,levels=c('EPDS','STAI-S','STAI-T'));v$model<-factor(v$model,levels=c('Baseline model','Full model (+ PD)'));assign(x,v)}
l$label<-gsub('Calibration intercept','Intercept',l$label,fixed=TRUE)
l$label<-gsub('Calibration slope','Slope',l$label,fixed=TRUE)
g<-ggplot(d,aes(x=predicted_prob))+
 geom_abline(intercept=0,slope=1,linetype='dotted',linewidth=.8)+
 geom_line(aes(y=observed_prob,color='Apparent',linetype='Apparent'),linewidth=1)+
 geom_line(aes(y=corrected_observed_prob,color='Optimism-corrected',linetype='Optimism-corrected'),linewidth=1.2)+
 facet_grid(outcome~model)+
 geom_text(data=l,aes(x=-Inf,y=Inf,label=label),inherit.aes=FALSE,hjust=-.025,vjust=1.15,size=4.1,lineheight=1.05)+
 coord_equal()+scale_linetype_manual(name=NULL,values=c(Apparent='dashed','Optimism-corrected'='solid'))+
 scale_color_manual(name=NULL,values=c(Apparent='#D55E00','Optimism-corrected'='#0072B2'))+
 labs(x='Predicted probability',y='Observed probability')+
 theme_classic(base_size=15)+theme(legend.position='bottom',legend.text=element_text(size=16),legend.key.width=grid::unit(1,'cm'),strip.text.x=element_text(size=14,face='bold'),strip.text.y=element_text(size=14,face='bold'),strip.background=element_rect(fill='white',colour='black'),axis.text=element_text(size=13),axis.title=element_text(size=15),plot.margin=margin(6,6,6,6))
ggsave(file.path(b,'Figure_2_larger_text.png'),g,width=8,height=10,dpi=400)
ggsave(file.path(b,'Figure_2_larger_text.tif'),g,width=8,height=10,dpi=400,compression='lzw')
