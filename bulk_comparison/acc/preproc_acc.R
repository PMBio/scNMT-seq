library(data.table)
library(purrr)
# library(rtracklayer)
library(Rsamtools)

############################
## Define I/O and options ##
############################

io <- list()
io$anno <- "/Users/ricard/data/NMT-seq/features/filt"
io$acc.data <- "/Users/ricard/data/encode/DNAse-seq/original/129_DHS_bamfile.bam"

opts <- list()
opts$chr <- c(1:19,"X","Y")

######################
## Load annotations ##
######################

anno = list.files(io$anno, pattern = ".bed") %>%
  lapply(function(x) {
    paste(io$anno, x,sep="/") %>%
      fread(select = c(1,2,3,5,6)) %>%
      setnames(c("chr", "start", "end","id","anno")) %>%
      .[, chr := gsub("chr", "", chr)]
  }) %>% rbindlist()


#########################
## Load DNase-seq data ##
#########################

# Define which options to extract from the BAM file
param = ScanBamParam(what = c("rname", "pos", "qwidth"))

# Parse BAM file
dhs_data = scanBam(io$acc.data, param = param) %>% .[[1]] %>% as.data.table() %>%
  setnames(c("chr", "start", "end")) %>%
  .[, chr := gsub("chr", "", chr)] %>% 
  .[, end := start + end] %>%
  .[chr %in% opts$chr] %>%
  .[complete.cases(.)]

total_reads = dhs_data[, .N]

# Overlap with annotations and calculate log2 
# (Q) Is there any relationship between number of reads and sequence properties? We should normalise by this
# (Q) Divide total number of reads as a library size correction? Why do we need to do this here?
# (Q) if you don't get any read, means the region was closed it was not captured?
# (Q) What about simple stochasticity in the number of reads? You could have an open region and get very few reads from it
dhs_data2 = setkey(dhs_data) %>%
  foverlaps(anno %>% setkey(chr, start, end), nomatch = 0) %>%
  # .[, .(dnase = .N *1e6 / total_reads), .(chr, start, end, anno, id)] %>%
  .[, .(dnase = .N), .(chr, start, end, anno, id)]

# (Q) Divide by length ???
dhs_data3 <- dhs_data2[,dnase_norm:=dnase/abs(end-start)] %>% 
  .[dnase_norm<50] %>%
  .[, dnase_log := log2(dnase_norm)]

# (Q) What is the differenc ebetween Filtered vs unfiltered alignment (.bam) files in the encode webpage???

# fwrite(dhs_data3, "/Users/ricard/data/encode/parsed.txt", quote = F, sep = "\t", row.names = F, col.names = T)
