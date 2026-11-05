# Called only after a successful correlation rerun; never reapproves model checkpoints.
.manifest_path <- here::here('output','corrected_file_manifest.csv')
.m <- utils::read.csv(.manifest_path,stringsAsFactors=FALSE)
.is_correlation <- .m$role=='correlation'
.m$md5[.is_correlation] <- unname(tools::md5sum(file.path(here::here(),.m$relative_path[.is_correlation])))
utils::write.csv(.m,.manifest_path,row.names=FALSE)
rm(.manifest_path,.m,.is_correlation)
