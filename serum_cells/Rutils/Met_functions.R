loadMultipleMetData <- function(files, context, osx=T) {
  if (osx) {
    data <- files %>% map(~ fread(sprintf("zcat < %s",.x), sep="\t", verbose=F, stringsAsFactors=F, showProgress=F))
  } else {
    data <- files %>% map(~ fread(sprintf("zcat %s",.x), sep="\t", verbose=F, stringsAsFactors=F, showProgress=F))
  }
  names(data) <- opts$cells
  list(data,names(data)) %>% pwalk(~.x[,sample:=.y])
  data <- rbindlist(data)
  colnames(data) <- c("chr","pos","rate","sample")
  return(data)
}