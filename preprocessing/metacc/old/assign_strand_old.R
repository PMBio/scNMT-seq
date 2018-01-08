# Script to: 
# - assign strand information to CpG sites from bismark file
# - remove non-interesting dinucleotides (non-cg methylation for example)
# - remove non-interesting chromosomes

# Input: 
# single-cell methylation files output from Bismark. In either one of the two following formats:
# input_format=1:
# chr     pos     rate
# 1       3019021 0
# 1       3027398 100
# 1       3052955 100

# input_format=2:
# chr     pos     met_reads non nomet_reads
# 1       3019021 0
# 1       3027398 100
# 1       3052955 100

# Output:
# output_format=1: we include the strand information but we don't modify the coordinates
# chr     pos     rate	strand
# 1       3019021 0	+
# 1       3027398 100	-
# 1       3052955 100	-

# output_format=2: we modify the coordinates to map negative CpGs to the positive strand 
# chr     pos     rate
# 1       3019021 0
# 1       3027397 100
# 1       3052954 100

suppressMessages(library(BSgenome.Mmusculus.UCSC.mm10))
suppressMessages(library(doParallel))
suppressMessages(library(data.table))
suppressMessages(library(purrr))
suppressMessages(library(Biostrings))

# Define I/0
io <- list()
if (grepl("[k|K]vothe",Sys.info()['nodename'])) {
  io$indir <- "/Users/ricard/data/NMT-seq/met/raw/allele_unspecific/stranded"
  io$outdir <- "/Users/ricard/data/NMT-seq/met/raw/allele_unspecific/unstranded"
} else {
  io$indir <- "/hps/nobackup/stegle/users/ricard/NMT-seq/met/raw/allele_unspecific/stranded"
  io$outdir <- "/hps/nobackup/stegle/users/ricard/NMT-seq/met/raw/allele_unspecific/unstranded"
}
dir.create(io$outdir)

# Define options
opts <- list()
opts$input_format <- 1
opts$output_format <- 2
opts$cores <- 3
opts$chr_list <- c(1:19,"X","Y")
# opts$remove_dinucleotides <- "non_cg"
opts$remove_dinucleotides <- NULL

# Load DNA files
# cat("Loading DNA files")
# dna.files <- list.files(io$indir.dnaseq,pattern="(.fa.gz)$")
# genome <- readDNAStringSet(paste(io$indir.dnaseq,dna.files,sep="/"), "fasta")
# names(genome) <- sapply(strsplit(names(genome), split=" "),"[[",1)

# Process samples
# samples <- sub(".tsv.gz","",list.files(io$indir,pattern="(.tsv.gz)$"))
# samples <- c("A10","B10","C10","D10","E10","F10","G10","H10","A11","B11","C11","D11","E11","F11","G11","H11","B12","C12","D12","E12","F12","G12")
samples <- c("E12","F12","G12")

registerDoParallel(cores=opts$cores)
invisible(foreach(i=1:length(samples)) %dopar% {
  
  # Load data
  print(samples[i])
  data <- fread(sprintf("zcat < %s/%s.tsv.gz",io$indir,samples[i]), sep="\t", header=T, verbose=F, showProgress=F)
  
  # Input format 1 (chr,pos,rate)
  if (opts$input_format == 1) {
    colnames(data) <- c("chr","pos","rate")
    
    # Input format 2 (chr,pos,met_reads,nonmet_reads)
  } else if (opts$input_format == 2) {
    colnames(data) <- c("chr","pos","met_reads","nonmet_reads")
    data[,rate:=round((met_reads/(met_reads+nonmet_reads))*100)]
  }
  
  # Remove weird chromosomes
  data <- data[chr %in% opts$chr_list,]
  
  # Get sequence
  seq <- getSeq(Mmusculus, paste0("chr",data$chr), data$pos-1, data$pos+1)
  data[,c("base_up","base","base_down") := list(substr(as.character(seq),1,1),substr(as.character(seq),2,2),substr(as.character(seq),3,3))]
  
  # Do sanity checks
  stopifnot(unique(data$base) %in% c("G","C"))
  
  # Remove dinucleotides (generally non-CpG methylation)
  if (length(opts$remove_dinucleotides)>0) {
    data[,dinucleotide:=ifelse(base=="C",paste0(base,base_down),paste0(base_up,base))]
    if (opts$remove_dinucleotides == "non_cg") {
      idx_keep <- which(data$dinucleotide=="CG")
      data <- data[idx_keep][,dinucleotide:=NULL]
      # Do sanity check
      data %>% split(.$base) %>% walk2(.,names(.), function(x,y) if (y=="C") { stopifnot(all(x$base_down=="G")) } else if (y=="G") { stopifnot(all(x$base_up=="C")) })
    } else {
      idx_keep <- which(!data$dinucleotide %in% opts$remove_dinucleotides)
      data <- data[idx_keep][,dinucleotide:=NULL]
    }
  }
  
  # Add strand information and remove base information
  data[,strand:=ifelse(base=="C","+","-")] %>% .[,c("base_up","base","base_down"):=NULL]
  
  # "Positivise" the dinucleotides and remove strand information
  if (opts$output_format == 2) {
    data[,pos:=ifelse(strand=="+",pos,pos-1)] %>% .[,strand:=NULL]
    # Sometimes we collect the same CpG site in both strands, in that case we can merge the information
    if (opts$input_format == 1) {
      data <- data[,.(rate=as.integer(round(mean(rate)))), keyby=c("chr","pos")]
    } else if (opts$input_format == 2) {
      data <- data[,.(rate=mean(rate), met_reads=sum(met_reads), nonmet_reads=sum(nonmet_reads)), keyby=c("chr","pos")]
    }
  }
  
  # Sanity check
  stopifnot(all(duplicated(data,by=c("chr","pos"))==FALSE))
  
  # Save results
  fwrite(data, file=sprintf("%s/%s.tsv",io$outdir,samples[i]), sep="\t", showProgress=FALSE, verbose=FALSE, col.names=FALSE)
})

# system(sprintf("gzip %s/*.tsv",io$outdir))
