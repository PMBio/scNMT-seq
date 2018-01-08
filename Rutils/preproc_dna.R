##################################################
##  First script to parse DNA accesibility data ##
##################################################

# This script performs quantification of features within samples

# It performs the following steps:
# - Manual sample filtering
# - Preprocessing of annotations: collect all GpC sites
# - Preprocessing of samples: collect all CpG sites from mm10 using the package "BSgenome.Mmusculus.UCSC.mm10"
# - Annotate samples with the preprocessed annotations

# Input:
# (1) two .gz compressed methylation files per sample with format  
# chr     pos      met_reads    nonmet_reads
# 19    3152031     1       0
# 19    3152424     0       1
# (2) genomic context annotation files in BED3 format
# (3) (opt) txt file with the samples to keep. If not provided, all samples are kept

# Output:
# - a tmp folder with a punch of rds files with the preprocessed samples and annotations. 
#   For example:
#     tmp/annos/genebody.rds
#     tmp/samples/3289STDY6312184.rds
#     tmp/annos_samples/3289STDY6312184_active_enhancers.rds
# - all.rds: all annotations and samples in one dataframe

options(warn=-1)
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(GenomicRanges))
suppressMessages(library(BSgenome.Mmusculus.UCSC.mm10))
suppressMessages(library(seqinr))
suppressMessages(library(stringr))
suppressMessages(library(doParallel))


#####################
## Define options ##
####################

## I/O ##

io <- list()

if (grepl("[k|K]vothe",Sys.info()['nodename'])) {
  # stop("Run this in the cluster")
  io$basedir <- "/Users/ricard/data/gastrulation/met"
} else {
  io$basedir <- "/hps/nobackup/stegle/users/ricard/gastrulation/met"
}


io$anno.folder <- "/Users/ricard/data/gastrulation/annotations/filt"
io$in.metadata <- str_c(io$anno.folder,"/metadata.txt")
io$in.folder <- str_c(io$basedir,"/unparsed")
io$out.folder <- str_c(io$basedir,"/tmp")
io$samples.file <- str_c(io$basedir,"/samples.txt")

## Options ##

opts <- list()

# Annotations to analyse
opts$anno <- "genebody"
# opts$anno <- "all"
# opts$anno <- c("3utr","5utr","cds_genebody","intergenic","noncds_genebody","prom2k", "prom2k_cgi", "prom2k_noncgi",
# "active_enhancers","primed_enhancers","poised_enhacners","LMRs","Tet1","Tet2","CGI","super_enhancers","p300","IAP")
if (opts$anno == "all") 
  opts$anno <- sapply(str_split(list.files(io$anno.folder, pattern = "\\.bed$"),"\\."),"[[", 1)

# Number of cores to use to parallelise the work
opts$cores <- 25

# Define valid chromosomes
opts$chr_list <- c(1:19,"X","Y")

cat("Processing methylation samples with the following options:\n")
cat(sprintf("- Input folder for annotation: %s\n",io$anno.folder))
cat(sprintf("- Input folder for methylation files: %s\n",io$in.folder))
cat(sprintf("- Output folder: %s\n",io$out.folder))
cat(sprintf("- Annotations: %s\n", paste(opts$anno, collapse=" ")))
cat(sprintf("- Number of cores: %d\n",opts$cores))
cat(sprintf("- Valid chromosomes: %s\n",str_c(opts$chr_list, collapse=" ")))
cat("\n")

###############
## Load data ##
###############

# Load samples to be kept
# TO-DO: ALL SAMPLES IF THIS NOT PROVIDED
samples_filt <- read.table(io$samples.file, header=F, stringsAsFactors=F)[,1]
stopifnot(all(!duplicated(samples_filt)))

# Load annotation metadata 
metadata <- read.table(io$in.metadata, header=T, stringsAsFactors=F)
# metadata <- read.table(io$in.metadata, header=T, stringsAsFactors=F)
colnames(metadata) <- c("chr","start","end","id","name")

############################
## Preprocess annotations ##
############################

cat("\nProcessing annotations...\n")

# Create ouput folder for annotations
dir.create(sprintf("%s/annos",io$out.folder), recursive=T)

# Run in parallel
source(str_c(io$basedir,"/unprocessed_data/test/utils.R"))
registerDoParallel(cores=opts$cores)
anno_list <- foreach(i=1:length(opts$anno)) %dopar% {
  anno=opts$anno[i]
  fname.out <- sprintf("%s/annos/%s.rds",io$out.folder,anno)
  if (file.exists(fname.out)) {
    cat(sprintf("annotation '%s' has already been processed, loading...\n",anno))
    granges <- readRDS(fname.out)
  } else {
    cat(sprintf("annotation '%s' has not been processed, processing...\n",anno))
    
    # Read annotation file
    anno.file <- sprintf("%s/%s.bed",io$anno.folder,anno)
    dat <- fread(anno.file ,sep="\t", header=F, select=c(1,2,3,4), verbose=F) %>% tbl_df
    colnames(dat) <- c("chr","start","end","id")
    stopifnot(all(dat$id %in% metadata$id))
    
    # Make sure there are no weird chromosomes
    dat$chr <- factor(convert_chr_format(dat$chr,to="short"))
    dat <- dat %>% filter(chr %in% opts$chr_list)
    
    # Collect all CpG sites for the corresponding annotation
    seq <- getSeq(Mmusculus, convert_chr_format(dat$chr,to="long"), dat$start,dat$end+1)
    matches <- vmatchPattern("CG",seq)
    CpG_sites <- start(matches) + dat$start - 1
    names(CpG_sites) <- dat$id
    
    # Count the number of CpG sites for each genomic context
    dat$n_CpG <- sapply(CpG_sites, length)
    
    # Create a GRanges object to do the overlap
    granges = makeGRangesFromDataFrame(dat, start.field="start", end.field="end",
                                       seqnames.field="chr", keep.extra.columns=TRUE)
    # Save the GRanges object
    saveRDS(granges,fname.out)
  }
  granges
}
names(anno_list) <- opts$anno


#########################################
## Preprocess and annotate the samples ##
#########################################

cat("\nAnnotating samples...\n")

# Read file names
filenames <- list.files(io$in.folder,pattern="(txt.gz)$")
stopifnot(all(samples_filt %in% unique(unname(sapply(filenames, function(f) strsplit(f,"_")[[1]][1])))))

# Create ouput folder for samples
dir.create(sprintf("%s/samples",io$out.folder), recursive=T)
dir.create(sprintf("%s/annos_samples",io$out.folder), recursive=T)

# Run in parallel
registerDoParallel(cores=opts$cores)
invisible(foreach(i=1:length(samples_filt)) %dopar% {
  sample=samples_filt[i]
  samples_processed <- list.files(sprintf("%s/samples",io$out.folder))
  if (sum(str_detect(samples_processed,str_c(sample,"_",opts$anno))) == length(opts$anno)) {
    cat(sprintf("Sample %s already processed for all required annotations...\n",sample)) 
  } else {
    fname.out <- sprintf("%s/samples/%s.rds",io$out.folder,sample)
    if (file.exists(fname.out)) {
      cat(sprintf("Sample %s has already been processed, loading...\n",sample))  
      dat <- readRDS(sprintf("%s/samples/%s.rds",io$out.folder,sample))
    } else {
      cat(sprintf("Sample %s has not been processed, loading...\n",sample))  
      
      # Read and parse raw methylation data
      fname.in <-  str_c(io$in.folder,"/",filenames[str_detect(filenames,sample)])
      dat1 <- fread(paste("gunzip -c", fname.in[1]), sep="\t", header=F, select=c(1,2,3,4), stringsAsFactors=TRUE, verbose=F) %>% tbl_df
      dat2 <- fread(paste("gunzip -c", fname.in[2]), sep="\t", header=F, select=c(1,2,3,4), stringsAsFactors=TRUE, verbose=F) %>% tbl_df
      colnames(dat1) <- c("chr","pos","met_reads","nonmet_reads")
      colnames(dat2) <- c("chr","pos","met_reads","nonmet_reads")
      
      # Merge the two files
      dat <- rbind(dat1,dat2) %>% arrange(chr,pos) %>% 
        group_by(chr,pos) %>% summarise(met_reads=sum(met_reads), nonmet_reads=sum(nonmet_reads)) %>% ungroup
      
      # Remove weird chromosomes
      dat <- dat %>% filter(chr %in% opts$chr_list)
      
      # Calculate methylation status for the CpG sites
      # dat$met_rate <- dat$met_reads/(dat$met_reads+dat$nonmet_reads)
      # dat$met_rate <- round(dat$met_rate,2)
      
      # Save the file for fast loading 
      saveRDS(dat,fname.out)
    }
    
    # Create GRanges objects for the overlap
    gr = makeGRangesFromDataFrame(dat, start.field="pos", end.field="pos", seqnames.field="chr", keep.extra.columns=TRUE)
    
    # Overlap data with annotations 
    for (anno in opts$anno) {
      if (file.exists(sprintf("%s/annos_samples/%s_%s.rds",io$out.folder,sample,anno))) {
        cat(sprintf("Annotation for %s with %s already found, skipping...\n",sample,anno)) 
      } else {
        cat(sprintf("Annotating %s with %s annotations...\n",sample,anno))
        
        # Overlap
        ov = findOverlaps(query=gr, subject=anno_list[[anno]], type="within", ignore.strand=TRUE)
        
        # Create dataframe with the methylation status of each independent CpG site
        shits <- anno_list[[anno]][subjectHits(ov),]
        qhits <- gr[queryHits(ov),]
        tmp <- data.frame(id=anno_list[[anno]]$id, n_CpG=anno_list[[anno]]$n_CpG)
        ov.df <- data.frame(id=shits$id, chr=seqnames(shits),
                            # start=start(ranges(shits)), end=end(ranges(shits)),
                            CpG_pos=start(ranges(qhits)),
                            met_reads=qhits$met_reads, nonmet_reads=qhits$nonmet_reads) %>% tbl_df %>% left_join(tmp, by="id")
        
        # Calculate methylation status for each region in the annotation by summarising over all CpG sites
        out <- ov.df %>% group_by(id) %>% 
          summarise(sample=sample, name=anno, 
                    chr=unique(chr), #start=unique(start), end=unique(end),
                    met_reads=sum(met_reads), nonmet_reads=sum(nonmet_reads), 
                    met_rate=round(met_reads/(met_reads+nonmet_reads),2),
                    weight=met_reads+nonmet_reads) %>% arrange(chr) # arrange(chr,start,end)
        
        # Save results
        saveRDS(out,sprintf("%s/annos_samples/%s_%s.rds",io$out.folder,sample,anno))
      }
    }
  }
})

# Join all annotated samples into one dataframe and save it
files <- list.files(str_c(io$out.folder,"/annos_samples"), pattern=".rds")
foo <- rbindlist(lapply(files, function(f) data.table(readRDS(str_c(io$out.folder,"/annos_samples/",f)))))
saveRDS(foo,str_c(io$out.folder,"/all.rds"))
