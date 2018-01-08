#############################################
## Second script to parse methylation data ##
#############################################

# This script performs basic processing of the methylation data:
# - Filter sites:
#     by minimum amount of reads in a minimum fraction of cells
# - Filter samples 
#     outliers
#     minimum coverage
# - (Opt) Merge sites with the same ens_id. For example, all exons in one gene.
#   The problem of this is that it looses the meaning of the start and end coordinates


options(warn=-1)
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(stringr))
suppressMessages(library(weights))
suppressMessages(library(argparse))

# Initialize argument parser
p <- ArgumentParser(description='')
# Define the fraction of sites to keep based on variance (1.0 to keep all sites)
p$add_argument('-c','--context', type="character",help='cg/CG or gc/GC')
p$add_argument('-a','--allele', action="store_true",help='Do allele-specific analysis?')
p$add_argument('-i','--imputed', action="store_true",help='Use imputed data?')

# Read arguments
args <- p$parse_args(commandArgs(TRUE))

####################
## Define options ##
####################

io <- list()
opts <- args

# Define what context to look at: CG (MET) or GC (DNA)
opts$context <- toupper(opts$context)
stopifnot(opts$context %in% c("CG","GC"))
stopifnot(opts$allele %in% c(FALSE,TRUE))

## I/O ##

# Own computer
if (grepl("[k|K]vothe",Sys.info()['nodename'])) {
  stop()
  # io$in.sample_metadata <- "/Users/ricard/NOMe-seq/data/sample_info.txt"
  # io$in.data <- "/Users/ricard/NOMe-seq/data/met/allele_inspecific/all.rds"
  # io$out.dir <- "/Users/ricard/NOMe-seq/data/met/allele_inspecific/filt"
  
# Cluster
} else {
  io$basedir <- "/hps/nobackup/stegle/users/ricard/NMT-seq"
  io$in.sample_metadata <- "/hps/nobackup/stegle/users/ricard/NMT-seq/sample_info.txt"
  if (opts$context == "GC") {
    if (opts$imputed) { cat("No imputed data for GC files"); stop() }
    if (opts$allele) {
      stop()
      io$in.data <- str_c(io$basedir,"/dna/parsed/allele_specific/unfilt/all.rds")
      io$out.dir <- str_c(io$basedir,"/dna/parsed/allele_specific/filt")
    } else {
      io$in.data <- str_c(io$basedir,"/dna/parsed/allele_unspecific/joint_module/unfilt/all.rds")
      io$out.dir <- str_c(io$basedir,"/dna/parsed/allele_unspecific/joint_module/filt")
    }
  } else {
    if (opts$allele) {
      if (opts$imputed) cat("No imputed data for allele_specific files"); stop()
      io$in.data <- str_c(io$basedir,"/met/parsed/allele_specific/unfilt/all.rds")
      io$out.dir <- str_c(io$basedir,"/met/parsed/allele_specific/filt")
    } else {
      if (opts$imputed) {
        io$in.data <- str_c(io$basedir,"/met/parsed/allele_unspecific/imputed/joint_module/unfilt/all.rds")
        io$out.dir <- str_c(io$basedir,"/met/parsed/allele_unspecific/imputed/joint_module/filt")
      } else {
        io$in.data <- str_c(io$basedir,"/met/parsed/allele_unspecific/unfilt/all.rds")
        io$out.dir <- str_c(io$basedir,"/met/parsed/allele_unspecific/filt")
      }
    }
  }
}


## Processing options ##

opts$filter <- F
opts$min_reads <- 1
opts$feature.min_cov <- 0.5

# Merge sites that have the same ens_id? (only available for some annotations)
opts$merge <- F
# opts$merge_annos <- c("5utr","3utr","cds_genebody","noncds_genebody","prom","prom_cgi","prom_noncgi","prom_active","prom_inactive","prom_poised")

# Annotations to analyse
# opts$annos <- c("CGI","Nanog","active_enhancers","primed_enhancers","prom_inactive","CTCF","Oct4","genebody","prom","prom_noncgi",
                # "DHS","Tet1","intergenic","prom_active","prom_poised","IAP","Tet2","p300","prom_cgi","super_enhancers")
opts$annos <- "all"


# Load samples to be kept
if (opts$context=="CG") {
  samples_keep <- read.table(io$in.sample_metadata, header=T, stringsAsFactors=F) %>% filter(passQC_met==T) %>% .$sample
} else if (opts$context=="GC") {
  samples_keep <- read.table(io$in.sample_metadata, header=T, stringsAsFactors=F) %>% filter(passQC_dna==T) %>% .$sample
}
if (opts$allele)
  samples_keep <- unlist(lapply(samples_keep, str_c, ".",c("genome1","genome2")))

cat("\nProcessing methylation samples with the following options:\n")
cat(sprintf("- Input data: %s\n",io$in.data))
cat(sprintf("- Output folder: %s\n",io$out.dir))
cat(sprintf("- Annotations: %s\n", paste(opts$annos, collapse=" ")))
cat(sprintf("- Allele-specific analysis: %s\n", as.character(opts$allele)))
cat(sprintf("- Annotating CG or GC?: %s\n",opts$context))
cat(sprintf("- Using CG imputed data: %s\n", as.character(opts$imputed)))
cat(sprintf("- Processing samples: %s\n",str_c(samples_keep, collapse=" ")))
cat("\n")

##################
## Loading data ##
##################

cat("Loading data...")

# Load data
data <- readRDS(io$in.data) %>% tbl_df

# Filter samples
data <- data %>% filter(sample %in% samples_keep)

# Filter annotations
if (opts$annos == "all") {
  opts$annos <- unique(data$anno)
} else {
  data <- filter(data, anno %in% opts$annos)
}

###################
## Sanity checks ##
###################

cat("Performing sanity checks...\n")

# Check that all annotations contain the same number of samples
tmp <- data %>% group_by(anno) %>% summarise(nsamples=length(unique(sample)))
stopifnot(length(unique(tmp$nsamples))==1)

#####################################
## Merge site with the same ens_id ##
#####################################

if (opts$merge) {
  cat("Merging sites with the same ensembl id...\n")
  
  data1 <- filter(data, anno %in% opts$merge_annos)
  data2 <- filter(data, !anno %in% opts$merge_annos)
  
  data1$ens_id <- unlist(lapply(strsplit(as.character(data1$id),"_"),"[[",1))
  
  data1 <- data1 %>% 
    group_by(chr,sample,anno,ens_id) %>% 
    summarise(rate=wtd.mean(rate,weight), weight=sum(weight)) %>% 
    ungroup %>% 
    dplyr::rename(id=ens_id)
  
  data <- rbind(data1,data2) %>% arrange(chr,id,anno)
}

################
## Parse data ##
################

# Calculate M value of methylation
offset <- 0.01
data <- data %>% mutate(m=log2((offset+rate/100)/(offset + 1-rate/100)))

###############
## Filtering ##
###############

if (opts$filter) {
  cat("Filtering...\n")
  
  ## Sites ##
  
  # Filter sites by minimum amount of reads opts$min_reads in a given fraction of cells opts$min_cov
  ncells <- length(unique(data$sample))
  tmp <- data %>% group_by(anno,id) %>% summarise(cov=sum(weight>=opts$min_reads)/ncells)
  filt_id <- tmp %>% filter(cov >= opts$feature.min_cov) %>% ungroup %>% .$id
  data <- data %>% filter(id %in% filt_id)
  
  # print results
  results <- tmp %>% group_by(anno) %>% summarise(good=sum(cov>=opts$feature.min_cov), bad=sum(cov<opts$min_cov), total=n())
  for (i in 1:nrow(results)) {
    anno <- results$anno[i]
    good <- results$good[i]
    total <- results$total[i]
    cat(sprintf("%s: %d/%d\n",anno,good,total))
  }
} else {
  cat("No filtering performed...\n")
}

## Cells ##

# Filter cells with a coverage less than opts$cell_mincov
# tmp <- data %>% group_by(sample) %>% summarise(coverage=n()/nrow(anno_metadata)) %>% arrange(coverage)
# hist(tmp$coverage, xlab="Coverage", main="Cell coverage")
# samples_filt <- tmp[tmp$coverage > opts$cell.min_cov,] %>% .$sample
# data <- data %>% filter(sample %in% samples_filt)

##################
## Save results ##
##################

cat("Saving results...\n")
dir.create(io$out.dir,showWarnings=F)

if (opts$context=="CG") {
  saveRDS(data,sprintf("%s/data.rds",io$out.dir))
} else if (opts$context=="GC") {
  saveRDS(data,sprintf("%s/data.rds",io$out.dir))
}

