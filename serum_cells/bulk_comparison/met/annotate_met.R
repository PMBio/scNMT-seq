###################################################################
##  Script to overlap bismark output files with genomic features ##
###################################################################

suppressMessages(library(data.table))
suppressMessages(library(seqinr))
suppressMessages(library(stringr))

## I/O ##
io <- list()
io$anno.folder <- "/Users/ricard/data/NMT-seq/features/filt"
io$in.folder <- "/Users/ricard/data/Ficz_2013/BS/raw"
io$out.folder <- "/Users/ricard/data/Ficz_2013/BS/parsed"


## Options ##
opts <- list()
opts$chr_list <- c(1:19,"X","Y") # Valid chromosomes
opts$annos <- "all" # Annotations to analyse
opts$samples <- c("rep1_parsed_mm10","rep2_parsed_mm10","rep3_parsed_mm10") # Samples

############################
## Preprocess annotations ##
############################

if (opts$annos == "all")
  opts$annos <- sapply(str_split(list.files(io$anno.folder, pattern = "\\.bed$"),"\\.bed"),"[[", 1)

# Run in parallel
anno_list <- list()
for (i in 1:length(opts$annos)) {
  
  # Read annotation file
  anno.file <- sprintf("%s/%s.bed",io$anno.folder,opts$anno[i])
  dat_anno <- fread(anno.file ,sep="\t", header=F, select=c(1,2,3,4,5), verbose=F) %>% 
    setnames(c("chr","start","end","strand","id"))
  
  # Check that there are no weird chromosomes
  anno_list[[i]] <- dat_anno %>% .[chr%in%opts$chr_list,] %>% setkey(chr,start,end)
}
names(anno_list) <- opts$anno



#########################################
## Preprocess and annotate the samples ##
#########################################

# Create ouput temporary folder
dir.create(sprintf("%s/tmp",io$out.folder), recursive=T)

# Run in parallel
for (i in 1:length(opts$samples)) {
  sample <- opts$samples[i]

  # Read and parse raw methylation data
  dat_sample <- fread(sprintf("zcat < %s/%s.txt.gz",io$in.folder,sample), sep="\t", verbose=F, showProgress=F) %>%
    .[,c("chr","pos","rate")]
  
  # Add 'start' and 'end' columns to do the overlap
  dat_sample <- dat_sample[,c("start","end") := list(pos,pos)][,pos:=NULL] %>% .[,chr:=as.factor(chr)] %>% setkey(chr,start,end)
  
  # Overlap data with annotations
  for (anno in opts$anno) {
    fname.out <- sprintf("%s/tmp/%s_%s.tsv",io$out.folder,sample,anno)
    if (file.exists(paste0(fname.out,".gz"))) {
      cat(sprintf("Annotation for %s with %s already found, loading...\n",sample,anno))
    } else {
      cat(sprintf("Annotating %s with %s annotations...\n",sample,anno))
      
      # Overlap
      ov <- foverlaps(dat_sample, anno_list[[anno]], nomatch=0) %>% .[,"i.end":=NULL] %>% setnames("i.start","pos")
      
      # Calculate methylation status for each region in the annotation by summarising over all CG sites
      # out <- ov[,c("sample","anno") := list(sample,anno)] %>% .[,.(rate=round(mean(rate)), var=round(var(rate)), weight=.N), keyby=.(sample,id,anno)]
      out <- ov[,c("sample","anno") := list(sample,anno)] %>% .[,.(rate=round(mean(rate)), weight=.N), keyby=.(sample,id,anno,chr,start,end)]
      
      # Store and save results
      fwrite(out, fname.out, quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)
    }
  }
}#)


# Compress output files
cat("Compressing...\n")
# system(sprintf("gzip -f %s/tmp/*.tsv",io$out.folder))
system(sprintf("pigz -f %s/tmp/*.tsv",io$out.folder))

# Concatenate everything and save it
for (i in opts$anno) {
  cat(i)
  files <- list.files(str_c(io$out.folder,"/tmp"), pattern=sprintf(".*%s.tsv.gz",i))
  foo <- lapply(files, function(f) fread(sprintf("zcat < %s/tmp/%s",io$out.folder,f))) %>% rbindlist
  write.table(foo, sprintf("%s/met_data_%s.tsv",io$out.folder,i), quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)
}

files <- list.files(str_c(io$out.folder,"/tmp"), pattern=".tsv.gz")
foo <- lapply(files, function(f) fread(sprintf("zcat < %s/tmp/%s",io$out.folder,f))) %>% rbindlist %>% setkey(anno)
write.table(foo, sprintf("%s/met_data.tsv",io$out.folder), quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)

system(sprintf("gzip -f %s/*.tsv",io$out.folder))
