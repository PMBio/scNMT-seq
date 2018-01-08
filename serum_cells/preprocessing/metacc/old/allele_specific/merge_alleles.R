
#######################################################
## script to merge allele-specific methylation files ##
#######################################################

suppressMessages(library(data.table))
suppressMessages(library(purrr))
suppressMessages(library(doParallel))
suppressMessages(library(argparse))

# Initialize argument parser
p <- ArgumentParser(description='')
p$add_argument('-c','--context', type="character",help='cg/CG or gc/GC')
p$add_argument('-n','--cores', type="integer",help='Number of cores')
opts <- p$parse_args(commandArgs(TRUE))

if (opts$context == "CG") {
  indir <- "/hps/nobackup/stegle/users/ricard/NMT-seq/met/raw/allele_specific"
  outdir <- "/hps/nobackup/stegle/users/ricard/NMT-seq/met/raw/allele_unspecific"
  samples <- unique(sapply(strsplit(list.files(indir,pattern="(csv.gz)$"),"_"), "[[",1))
} else if (opts$context == "GC") {
  indir <- "/hps/nobackup/stegle/users/ricard/NMT-seq/dna/raw/allele_specific"
  outdir <- "/hps/nobackup/stegle/users/ricard/NMT-seq/dna/raw/allele_unspecific"
  samples <- unique(sapply(strsplit(list.files(indir,pattern="(csv.gz)$"),"_"), "[[",1))
} else {
  stop("wrong context")
}
dir.create(outdir, showWarnings=FALSE)

registerDoParallel(cores=opts$cores)
invisible(foreach(i=1:length(samples)) %dopar% {
  sample=samples[i]
  outfile <- sprintf("%s/%s.tsv",outdir,sample)
  if (file.exists(outfile)) {
    print(sprintf("%s already processed...", sample))
  } else {
    print(sprintf("processing %s...", sample))
    genome1 <- fread(sprintf("zcat < %s/%s_genome1.csv.gz",indir,sample), sep=",", header=T, showProgress=F)[,c(1,2,5,6)]
    genome2 <- fread(sprintf("zcat < %s/%s_genome2.csv.gz",indir,sample), sep=",", header=T, showProgress=F)[,c(1,2,5,6)]
    unassigned <- fread(sprintf("zcat < %s/%s_unassigned.csv.gz",indir,sample), showProgress=F)[,c(1,2,5,6)]
    colnames(genome1) <- colnames(genome2) <- colnames(unassigned) <- c("chr","pos","met_reads","nonmet_reads")
    
    merged <- rbind(genome1,genome2,unassigned) %>% .[,list(met_reads=sum(met_reads), nonmet_reads=sum(nonmet_reads)), by=c("chr","pos")]
    
    fwrite(merged, outfile, sep="\t", showProgress=FALSE, verbose=FALSE, col.names=TRUE)
  }
})
system(sprintf("gzip -f %s/*.tsv",outdir))