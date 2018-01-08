###################################################################
##  Script to overlap bismark output files with genomic features ##
###################################################################

# This script overlaps the output bismark files (individual CpG sites) with genomic features such as promoters, gene bodies, etc.

# - Preprocessing of annotations: collect all CpG sites
# - Preprocessing of samples: collect all CpG sites from mm10 using the package "BSgenome.Mmusculus.UCSC.mm10"
# - Annotate samples with the preprocessed annotations

## Input ##
# (1) output bismark file with one of the two following formats
# input_format=1
#   chr     pos      met_reads    nonmet_reads
#   19    3152031     1       0
#   19    3152424     0       1
# input_format=2
#   chr     pos      rate
#   19    3152031     100
#   19    3152424     0

# (2) genomic feature annotation files in BED6 format
#   1	3531624	3531843	*	CGI_1	CGI
#   1	3670619	3671074	*	CGI_2	CGI
#   1	3671654	3672156	*	CGI_3	CGI

## Output ##
# (1) a tmp folder with a punch of tsv files with the preprocessed samples and annotations. 
#   For example:
#     tmp/[SAMPLE]_[FEATURE].tsv
# (2) all.tsv: all annotations and samples in one dataframe
# output format=1
#   sample	id	anno	rate	weight
#   3289STDY6312493	super_enhancers_100	super_enhancers	42	31
#   3289STDY6312493	super_enhancers_1001	super_enhancers	0	1
#   3289STDY6312493	super_enhancers_1002	super_enhancers	0	2

# """
# To-do:
# - Currently only input_format=2 is accepted
# - Currently only output_format=1 is accepted
# - 
# """

options(warn=-1)
suppressMessages(library(data.table))
# suppressMessages(library(BSgenome.Mmusculus.UCSC.mm10))
suppressMessages(library(seqinr))
suppressMessages(library(stringr))
suppressMessages(library(doParallel))
suppressMessages(library(argparse))

# Initialize argument parser
p <- ArgumentParser(description='')
p$add_argument('-c','--context', type="character",help='cg/CG or gc/GC')
p$add_argument('-n','--cores', type="integer",help='Number of cores')

# Read arguments
args <- p$parse_args(commandArgs(TRUE))

#####################
## Define options ##
####################

# args <- list()
# args$context <- "GC"
# args$cores <- 4

io <- list()
opts <- args

## I/O ##

# Define what context to look at: CG (MET) or GC (ACC)
opts$context <- toupper(opts$context)
stopifnot(opts$context %in% c("CG","GC"))

## Own computer ##
if (grepl("ricard",Sys.info()['nodename'])) {
  io$basedir <- "/Users/ricard/data/NMT-seq_EB"
  io$anno.folder <- str_c(io$basedir,"/features/filt")
  io$in.sample_metadata <- str_c(io$basedir,"/sample_sheet.csv")
  
  # GC
  if (opts$context == "GC") {
    io$in.folder <- str_c(io$basedir,"/acc/raw/merged/binarised")
    io$out.folder <- str_c(io$basedir,"/acc/parsed")
  # CG
  } else if (opts$context=="CG") {
    io$in.folder <- str_c(io$basedir,"/met/raw/merged/binarised")
    io$out.folder <- str_c(io$basedir,"/met/parsed")
  }
  
## Cluster ##
} else {
  stop()
  io$basedir <- "/hps/nobackup/stegle/users/ricard/NMT-seq_EB"
  io$anno.folder <- str_c(io$basedir,"/features/filt")
  io$in.sample_metadata <- str_c(io$basedir,"/sample_info_all.txt")
  
  # GC
  if (opts$context == "GC") {
    io$in.folder <- str_c(io$basedir,"/acc/raw/filtered/unstranded/binarised")
    io$out.folder <- str_c(io$basedir,"/acc/parsed")
  # CG
  } else if (opts$context=="CG") {
    opts$sep <- "\t"
    io$in.folder <- str_c(io$basedir,"/met/raw/filtered/unstranded/binarised")
    io$out.folder <- str_c(io$basedir,"/met/parsed")
  }
}


## Options ##

# Valid chromosomes
opts$chr_list <- c(1:19,"X","Y")

# Annotations to analyse
opts$annos <- "all"
if (opts$annos == "all")
  opts$annos <- sapply(str_split(list.files(io$anno.folder, pattern = "\\.bed$"),"\\."),"[[", 1)

cat("\nProcessing methylation samples with the following options:\n")
cat(sprintf("- Input folder for annotation: %s\n",io$anno.folder))
cat(sprintf("- Input folder for bismark files: %s\n",io$in.folder))
cat(sprintf("- Output folder: %s\n",io$out.folder))
cat(sprintf("- Annotations: %s\n", paste(opts$annos, collapse=" ")))
cat(sprintf("- Number of cores: %d\n",opts$cores))
cat(sprintf("- Valid chromosomes: %s\n",str_c(opts$chr_list, collapse=" ")))
cat(sprintf("- Annotating CG or GC?: %s\n",opts$context))
cat("\n")

###############
## Load data ##
###############

# Load samples to be kept
if (opts$context=="CG") {
  samples_keep <- fread(io$in.sample_metadata, header=T) %>% .[pass_metQC==T,sample]
} else if (opts$context=="GC") {
  samples_keep <- fread(io$in.sample_metadata, header=T) %>% .[pass_accQC==T,sample]
} else{
  stop()
}
stopifnot(all(!duplicated(samples_keep)))

cat(sprintf("- Processing samples: %s\n",str_c(samples_keep, collapse=" ")))


############################
## Preprocess annotations ##
############################

# Run in parallel
# registerDoParallel(cores=args$cores)
# anno_list <- foreach(i=1:length(opts$anno)) %dopar% {
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
# registerDoParallel(cores=args$cores)
# invisible(foreach(i=1:length(samples_keep)) %dopar% {
for (i in 1:length(samples_keep)) {
  sample=samples_keep[i]
  samples_processed <- list.files(sprintf("%s/tmp",io$out.folder))
  if (sum(str_detect(samples_processed,str_c(sample,"_",opts$anno))) == length(opts$anno)) {
    cat(sprintf("Sample %s already processed for all required annotations...\n",sample)) 
  } else {
    cat(sprintf("Sample %s has not been processed, annotating...\n",sample))  
    
    # Read and parse raw methylation data
    dat_sample <- fread(sprintf("zcat < %s/%s.tsv.gz",io$in.folder,sample), sep="\t", header="auto", verbose=F, showProgress=F)
    colnames(dat_sample) <- c("chr","pos","rate")
    
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
        out <- ov[,c("sample","anno") := list(sample,anno)] %>% .[,.(rate=round(mean(rate)), var=round(var(rate)), weight=.N), keyby=.(sample,id,anno)]
        
        # Store and save results
        fwrite(out, fname.out, quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)
      }
    }
  }
}#)


# Compress output files
cat("Compressing...\n")
# system(sprintf("gzip -f %s/tmp/*.tsv",io$out.folder))
system(sprintf("pigz -p %d -f %s/tmp/*.tsv",args$cores, io$out.folder))

# Concatenate everything and save it
cat("Annotations finished, combining results...\n")
files <- list.files(str_c(io$out.folder,"/tmp"), pattern=".tsv.gz")
foo <- lapply(files, function(f) fread(sprintf("zcat < %s/tmp/%s",io$out.folder,f))) %>% rbindlist

if (opts$context == "CG") {
  write.table(foo, sprintf("%s/met_data.tsv",io$out.folder), quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)
} else {
  write.table(foo, sprintf("%s/acc_data.tsv",io$out.folder), quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)
}
system(sprintf("gzip -f %s/*.tsv",io$out.folder))
