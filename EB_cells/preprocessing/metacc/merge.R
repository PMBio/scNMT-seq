
##########################################
## Script to merge bismark output files ##
##########################################

# The output of bismark comes as two files per sample: R1 and R2
# This script merges them and generates a single file per sample.
# If the same site is present in the two files, then the reads are simply summed

library(data.table)
library(purrr)
library(doParallel)

indir <- "/Users/ricard/data/NMT-seq_EB/acc/raw"
outdir <- "/Users/ricard/data/NMT-seq_EB/acc/raw/merged"; dir.create(outdir)
cores <- 3

filenames <- list.files(indir,pattern="(tsv.gz)$")
# samples <- unique(substr(filenames,1,nchar(filenames)-10))
# samples <- unique(sapply(strsplit(filenames,split="_"),"[[",1))
samples <- unique( paste( sapply(strsplit(filenames,split="_"),"[[",1), sapply(strsplit(filenames,split="_"),"[[",2), sep="_" ) )


# registerDoParallel(cores=cores)
# invisible(foreach(i=1:length(samples)) %dopar% {
for (i in 1:length(samples)) {
  print(sprintf("%s (%d/%d)", samples[i], i, length(samples)))
  fname.in <- sprintf("%s/%s",indir,filenames[grep(samples[i],filenames)])
  
  if (length(fname.in) == 2) {
    dat1 <- fread(sprintf("zmore %s",fname.in[1]), header=F, select=c(1,2,3,4), stringsAsFactors=TRUE, verbose=F)
    dat2 <- fread(sprintf("zmore %s",fname.in[2]), header=F, select=c(1,2,3,4), stringsAsFactors=TRUE, verbose=F)
    colnames(dat1) <- c("chr","pos","met_reads","nonmet_reads")
    colnames(dat2) <- c("chr","pos","met_reads","nonmet_reads")
    dat <- rbind(dat1,dat2) %>% .[,.(met_reads=sum(met_reads), nonmet_reads=sum(nonmet_reads)), by=c("chr","pos")] %>% setkey(chr,pos)
    
    # In some cases, samples have four files (WHY?)
  } else if (length(fname.in) == 4) {
    dat1 <- fread(sprintf("zmore %s",fname.in[1]), header=F, select=c(1,2,5,6), stringsAsFactors=TRUE, verbose=F)
    dat2 <- fread(sprintf("zmore %s",fname.in[2]), header=F, select=c(1,2,5,6), stringsAsFactors=TRUE, verbose=F)
    dat3 <- fread(sprintf("zmore %s",fname.in[3]), header=F, select=c(1,2,5,6), stringsAsFactors=TRUE, verbose=F)
    dat4 <- fread(sprintf("zmore %s",fname.in[4]), header=F, select=c(1,2,5,6), stringsAsFactors=TRUE, verbose=F)
    colnames(dat1) <- c("chr","pos","met_reads","nonmet_reads")
    colnames(dat2) <- c("chr","pos","met_reads","nonmet_reads")
    colnames(dat3) <- c("chr","pos","met_reads","nonmet_reads")
    colnames(dat4) <- c("chr","pos","met_reads","nonmet_reads")
    dat <- rbind(dat1,dat2,dat3,dat4) %>% .[,.(met_reads=sum(met_reads), nonmet_reads=sum(nonmet_reads)), by=c("chr","pos")] %>% setkey(chr,pos)
    
  } else {
    stop("error")
  }
  fwrite(dat, file=sprintf("%s/%s.tsv",outdir,samples[i]), sep="\t", showProgress=FALSE, verbose=FALSE, col.names=TRUE)
}
# )

# system(sprintf("gzip -f %s/*.tsv",io$outdir))
system(sprintf("pigz -p %d -f %s/*.tsv", cores, outdir))