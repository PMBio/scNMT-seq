# Create a bulk genome

library(data.table)
library(purrr)

io <- list()
# io$indir <- "/hps/nobackup/stegle/users/ricard/NMT-seq/dna/raw/merged"
# io$outfile <- "/hps/nobackup/stegle/users/ricard/NMT-seq/dna/raw/merged/bulk.txt"
io$indir <- "/hps/nobackup/stegle/users/ricard/NMT-seq/met/raw/merged"
io$outfile <- "/hps/nobackup/stegle/users/ricard/NMT-seq/met/raw/merged/bulk.txt"

cat("Loading data...")
data <- paste(io$indir,list.files(io$indir),sep="/") %>% 
  map(~ fread(sprintf("zcat < %s",.x), sep="\t", header=T, verbose=F, stringsAsFactors=F, showProgress=F)) %>%
  rbindlist

cat("Processing...")
data_filt <- data[,list(rate=round(mean(rate)), weight=.N), by=.(chr,pos)]

cat("Saving...")
fwrite(data_filt, file = io$outfile, append = FALSE, quote = "auto",
       sep = "\t", row.names = FALSE, col.names = TRUE)
system(sprintf("gzip %s",io$outfile))
