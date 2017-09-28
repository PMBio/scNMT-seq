########################################################################
## Script to assign strand information to CpG sites from bismark file ##
########################################################################

# Input: 
# single-cell methylation files output from Bismark. In either one of the two following formats:
# input_format=1:
  # chr     pos     rate
  # 1       3019021 0
  # 1       3027398 100
  # 1       3052955 100

# input_format=2:
  # chr     pos     met_reads non nomet_reads
  # 1       3019021 0 1
  # 1       3027398 1 1
  # 1       3052955 1 0

# Output:
# output_format=1: we include the strand information in a new column but we don't modify the coordinates
# output_format=2: we modify the coordinates to map negative CpGs to the positive strand 

# Load libraries
suppressMessages(library(BSgenome.Mmusculus.UCSC.mm10))
suppressMessages(library(doParallel))
suppressMessages(library(data.table))
suppressMessages(library(purrr))
suppressMessages(library(Biostrings))
suppressMessages(library(argparse))

# Initialize argument parser
p <- ArgumentParser(description='')
p$add_argument('-c','--context', type="character",help='cg/CG or gc/GC')
p$add_argument('-n','--cores', type="integer",help='Number of cores')
opts <- p$parse_args(commandArgs(TRUE))
opts$context <- toupper(opts$context); stopifnot(opts$context %in% c("CG","GC"))

# Define options
opts$input_format <- 1
opts$output_format <- 2

# Define I/0
io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  stop()
} else {
  if (opts$context=="CG") {
    io$indir <- "/hps/nobackup/stegle/users/ricard/NMT-seq/met/raw/allele_unspecific/filtered"
    io$outdir <- "/hps/nobackup/stegle/users/ricard/NMT-seq/met/raw/allele_unspecific/filtered/unstranded"
  } else {
    io$indir <- "/hps/nobackup/stegle/users/ricard/NMT-seq/acc/raw/allele_unspecific/filtered"
    io$outdir <- "/hps/nobackup/stegle/users/ricard/NMT-seq/acc/raw/allele_unspecific/filtered/unstranded"
  }
}
dir.create(io$outdir)

# Load samples
samples <- sub(".tsv.gz","",list.files(io$indir,pattern="(.tsv.gz)$"))

# Parallelise processing
registerDoParallel(cores=opts$cores)
invisible(foreach(i=1:length(samples)) %dopar% {
  outfile <- sprintf("%s/%s.tsv",io$outdir,samples[i])
  if (file.exists(outfile)) {
    cat(sprintf("Sample %s already processed, skipping...\n",samples[i]))
  } else {
    cat(sprintf("Processing %s...\n",samples[i]))
    
    # Load data
    data <- fread(sprintf("zcat < %s/%s.tsv.gz",io$indir,samples[i]), verbose=F, showProgress=F)
    
    # Input format 1 (chr,pos,rate)
    if (opts$input_format == 1) {
      colnames(data) <- c("chr","pos","rate")
      
      # Input format 2 (chr,pos,met_reads,nonmet_reads)
    } else if (opts$input_format == 2) {
      colnames(data) <- c("chr","pos","met_reads","nonmet_reads")
      # data[,rate:=round((met_reads/(met_reads+nonmet_reads))*100)]
    }
    
    # Get sequence
    seq <- getSeq(Mmusculus, sub("MT","M",paste0("chr",data$chr)), data$pos-1, data$pos+1)
    data[,c("base_up","base","base_down") := list(substr(as.character(seq),1,1),substr(as.character(seq),2,2),substr(as.character(seq),3,3))]
    
    # Do sanity checks
    stopifnot(unique(data$base) %in% c("G","C"))
    
    # Add strand information and remove base information
    if (opts$context=="CG") {
      data[,strand:=ifelse(base=="C","+","-")] %>% .[,c("base_up","base","base_down"):=NULL]
    } else {
      data[,strand:=ifelse(base=="G","+","-")] %>% .[,c("base_up","base","base_down"):=NULL]
    }
    
    # "Positivise" the dinucleotides and remove strand information
    if (opts$output_format == 2) {
      data[,pos:=ifelse(strand=="+",pos,pos-1)] %>% .[,strand:=NULL]
      # Sometimes we collect the same CpG site in both strands, in that case we can merge the information
      if (opts$input_format == 1) {
        data <- data[,.(rate=as.integer(round(mean(rate)))), keyby=c("chr","pos")]
      } else if (opts$input_format == 2) {
        data <- data[,.(met_reads=sum(met_reads), nonmet_reads=sum(nonmet_reads)), keyby=c("chr","pos")]
      }
    }
    
    # Sanity check
    stopifnot(all(duplicated(data,by=c("chr","pos"))==FALSE))
    
    # Save results
    fwrite(data, outfile, sep="\t", showProgress=FALSE, verbose=FALSE, col.names=FALSE)
  }
})

# system(sprintf("gzip -f %s/*.tsv",io$outdir))
system(sprintf("pigz -p %d -f %s/*.tsv", opts$cores, io$outdir))
