# Reproduce corrected observed-data Tables 1a and 1b (28 October 2026).
# Called by section 2.2.1 of 3PAD_PLOS_revision.qmd after outcome construction.
# No model fitting, imputation or bootstrap execution occurs here.
rebuild_tables_1a_1b <- function(table1_data, output_dir) {
stopifnot(nrow(table1_data)==91L,
 all(c("country","pos_LE1","outcome_epds","outcome_stai_s","outcome_stai_t") %in% names(table1_data)),
 !anyNA(table1_data[c("outcome_epds","outcome_stai_s","outcome_stai_t")]),
 sum(table1_data$outcome_epds==1L)==13L,
 sum(table1_data$outcome_stai_s==1L)==12L,
 sum(table1_data$outcome_stai_t==1L)==21L)
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
d<-table1_data
vars<-c('country','gb1_age','edu','first_preg','bmi','workFT','pcos','exercise1','exercise2','SSRI','pos_LE1','pos_LE2','neg_LE1','neg_LE2','gb1_epds_sum','gb2_epds_sum','gb1_stai_s_sum','gb2_stai_s_sum','gb1_stai_t_sum','gb2_stai_t_sum')
res<-list();stats<-list()
for(outcome in c('outcome_epds','outcome_stai_s','outcome_stai_t'))for(v in vars){
 x<-d[[v]];g<-d[[outcome]];num<-is.numeric(x)
 for(k in 0:1){xx<-x[g==k];N<-length(xx);n<-sum(!is.na(xx));miss<-sum(is.na(xx))
  if(num)res[[length(res)+1L]]<-data.frame(outcome,variable=v,group=k,category='Mean (SD)',count=n,N,missing=miss,value=sprintf('%.1f (%.1f)',mean(xx,na.rm=T),sd(xx,na.rm=T)),mean=mean(xx,na.rm=T),sd=sd(xx,na.rm=T))
  else for(lv in levels(x))res[[length(res)+1L]]<-data.frame(outcome,variable=v,group=k,category=lv,count=sum(xx==lv,na.rm=T),N,missing=miss,value=sprintf('%d (%.1f)',sum(xx==lv,na.rm=T),100*sum(xx==lv,na.rm=T)/N),mean=NA,sd=NA)
 }
 if(num){
  sh<-shapiro.test(x)$p.value;welch<-t.test(x~g)$p.value;wil<-suppressWarnings(wilcox.test(x~g)$p.value)
  orig<-if(sh>.05)welch else wil;test<-if(sh>.05)'Welch t-test' else 'Wilcoxon rank-sum'
  reported<-welch
 }else{
  t<-table(x,g);fish<-fisher.test(t,workspace=2e8)$p.value;chi<-suppressWarnings(chisq.test(t)$p.value)
  orig<-if(any(t<5))fish else chi;test<-if(any(t<5))'Fisher exact' else 'Chi-square (Yates correction)';reported<-chi;sh<-NA;welch<-NA;wil<-NA
 }
 stats[[length(stats)+1L]]<-data.frame(outcome,variable=v,original_workflow_test=test,original_workflow_p=orig,manuscript_described_test_p=reported,shapiro_p=sh)
}
write.csv(do.call(rbind,res),file.path(output_dir,'table1_observed_summaries.csv'),row.names=F)
write.csv(do.call(rbind,stats),file.path(output_dir,'table1_test_audit.csv'),row.names=F)

summary_long <- do.call(rbind,res)
tests <- do.call(rbind,stats)
labels <- c(country="Place of birth",gb1_age="Maternal age",edu="Education",
 first_preg="First pregnancy",bmi="BMI",workFT="Employment or study status",
 pcos="PCOS",exercise1="Exercise early pregnancy",exercise2="Exercise late pregnancy",
 SSRI="SSRI use early pregnancy",pos_LE1="PLE early pregnancy",pos_LE2="PLE late pregnancy",
 neg_LE1="NLE early pregnancy",neg_LE2="NLE late pregnancy",
 gb1_epds_sum="EPDS early pregnancy",gb2_epds_sum="EPDS late pregnancy",
 gb1_stai_s_sum="STAI-S early pregnancy",gb2_stai_s_sum="STAI-S late pregnancy",
 gb1_stai_t_sum="STAI-T early pregnancy",gb2_stai_t_sum="STAI-T late pregnancy")
category_order <- list(country=c("Sweden","Other"),edu=c("University","Less"),
 first_preg=c("No","Yes"),workFT=c("Full-time","Other"),pcos=c("No","Yes"),
 exercise1=c("Low to none","Moderate to high"),exercise2=c("Low to none","Moderate to high"),
 SSRI=c("No","Yes"),pos_LE1=c("No","Yes"),pos_LE2=c("No","Yes"),
 neg_LE1=c("No","Yes"),neg_LE2=c("No","Yes"))
format_p <- function(p) if(p < .001) "<.001" else sub("^0", "", sprintf("%.3f",p))
get_summary <- function(o,v,k,cat) {
 z <- summary_long[summary_long$outcome==o & summary_long$variable==v &
  summary_long$group==k & summary_long$category==cat,,drop=FALSE]
 stopifnot(nrow(z)==1L);z
}
make_table <- function(outcomes) {
 rows <- list()
 value_names <- unlist(lapply(outcomes,function(o)paste0(o,c("_low","_high"))))
 p_names <- paste0(outcomes,"_p")
 for(v in vars) {
  cats <- category_order[[v]]
  if(is.null(cats))cats <- "Mean (SD)"
  for(cat in cats) {
   row <- list(Variable=unname(labels[[v]]),Category=if(cat=="Less")"No university" else cat)
   values <- unlist(lapply(outcomes,function(o)sapply(0:1,function(k) {
    z<-get_summary(o,v,k,cat)
    if(cat=="Mean (SD)" && z$missing>0)paste0(z$value,"; n = ",z$count) else z$value
   })),use.names=FALSE)
   row[value_names]<-as.list(values)
   row[p_names]<-lapply(outcomes,function(o) {
    z<-tests[tests$outcome==o & tests$variable==v,,drop=FALSE]
    if(cat==tail(cats,1))format_p(z$original_workflow_p) else ""
   })
   rows[[length(rows)+1L]]<-as.data.frame(row,check.names=FALSE,stringsAsFactors=FALSE)
  }
  missing_any<-any(summary_long$missing[summary_long$variable==v & summary_long$outcome %in% outcomes]>0)
  if(cats[1]!="Mean (SD)" && missing_any) {
   row<-list(Variable=unname(labels[[v]]),Category="Missing")
   row[value_names]<-as.list(unlist(lapply(outcomes,function(o)sapply(0:1,function(k) {
    z<-get_summary(o,v,k,cats[1]);sprintf("%d (%.1f)",z$missing,100*z$missing/z$N)
   })),use.names=FALSE))
   row[p_names]<-as.list(rep("",length(outcomes)))
   rows[[length(rows)+1L]]<-as.data.frame(row,check.names=FALSE,stringsAsFactors=FALSE)
  }
 }
 do.call(rbind,rows)
}
table1a <- make_table("outcome_epds")
table1b <- make_table(c("outcome_stai_s","outcome_stai_t"))
names(table1a) <- c("Variable","Category","EPDS <12 (n=78)","EPDS >=12 (n=13)","p")
names(table1b) <- c("Variable","Category","STAI-S <40 (n=79)","STAI-S >=40 (n=12)",
 "STAI-T <40 (n=70)","STAI-T >=40 (n=21)","STAI-S p","STAI-T p")
write.csv(table1a,file.path(output_dir,"Table_1a_corrected.csv"),row.names=FALSE,na="")
write.csv(table1b,file.path(output_dir,"Table_1b_corrected.csv"),row.names=FALSE,na="")
notes <- paste(
 "Values are n (%) or mean (SD). Percentages use the full outcome-group denominator, including missing entries.",
 "Categorical missing counts are explicit. Continuous summaries and tests use observed data, with available n shown when reduced.",
 "First pregnancy denotes pregnancy history, not parity. Employment/study includes full-time employment, self-employment or study.",
 "Tests retain the original workflow: Welch t-test if pooled Shapiro-Wilk p > .05, otherwise Wilcoxon rank-sum;",
 "Fisher exact if any observed categorical cell <5, otherwise chi-square with the default Yates correction.",
 "Tests exclude observations missing the relevant variable. Unadjusted p-values are rounded to three decimals; p < .001 is displayed as <.001.",
 sep="\n")
writeLines(notes,file.path(output_dir,"Table_1_notes.txt"))
if(requireNamespace("knitr",quietly=TRUE)) {
 html <- c("<!doctype html><html><head><meta charset='utf-8'><title>Corrected Tables 1a and 1b</title>",
 "<style>body{font-family:serif;margin:2em}table{border-collapse:collapse;margin-bottom:2em}th,td{padding:5px 10px;border-bottom:1px solid #ddd;text-align:left}</style></head><body>",
 "<h2>Table 1a. Participant characteristics by postpartum depressive symptoms</h2>",
 as.character(knitr::kable(table1a,format="html",escape=TRUE)),
 "<h2>Table 1b. Participant characteristics by postpartum state and trait anxiety</h2>",
 as.character(knitr::kable(table1b,format="html",escape=TRUE)),
 paste0("<p>",gsub("\n","<br>",notes),"</p></body></html>"))
 writeLines(html,file.path(output_dir,"Tables_1a_1b_corrected.html"))
}
list(table1a=table1a,table1b=table1b,summaries=summary_long,tests=tests,notes=notes)

}
