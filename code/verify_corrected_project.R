verify_corrected_project <- function(project_root=getwd()) {
  project_root<-normalizePath(project_root,mustWork=TRUE)
  manifest_path<-file.path(project_root,'output','corrected_file_manifest.csv')
  if(!file.exists(manifest_path))stop('Corrected-results manifest is missing. Do not use the archived results.')
  manifest<-utils::read.csv(manifest_path,stringsAsFactors=FALSE)
  paths<-file.path(project_root,manifest$relative_path)
  if(any(!file.exists(paths)))stop('A required corrected file is missing: ',manifest$relative_path[which(!file.exists(paths))[1]])
  actual<-unname(tools::md5sum(paths))
  if(!identical(actual,manifest$md5))stop('A corrected input/result file has changed. Verify or regenerate results before reusing checkpoints: ',manifest$relative_path[which(actual!=manifest$md5)[1]])
  for(key in c('epds','stais','stait')) {
    files<-list.files(file.path(project_root,'output/bootstrap',key),pattern=paste0('^',key,'_bootstrap_[0-9]{4}[.]rds$'))
    if(!identical(sort(files),sprintf(paste0(key,'_bootstrap_%04d.rds'),1:1000)))stop('The corrected ',key,' checkpoint set is incomplete.')
  }
  invisible(TRUE)
}
