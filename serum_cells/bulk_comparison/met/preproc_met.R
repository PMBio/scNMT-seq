library(data.table)
library(purrr)
library(rtracklayer)

io <- list()
io$met.data <- "/Users/ricard/data/Ficz_2013/BS/E14_ES_Serum_P12_3_CpG_methylation_calls.txt.gz"

# Load bulk bisulfie data 
bulk_bs = io$met.data %>%
  paste0("zcat < ",.) %>%
  fread(select = 2:4, fill = TRUE) %>%
  setnames(c("rate", "chr", "start")) %>%
  .[rate == "+", rate := "100"] %>%
  .[rate == "-", rate := "0"] %>%
  .[, rate := as.numeric(rate)] %>%
  .[!is.na(rate)] %>%
  .[,.(rate=mean(rate), N=.N),by=c("chr","start")] %>%
  .[,rate:=round(rate,2)] %>%
  setkey(chr,start)

# Liftover to mm10
# A chain file essentially details many local alignments, so it is
# possible for the "from" ranges to map to overlapping regions in
# the other sequence. The "from" ranges are guaranteed to be
# disjoint (but do not necessarily cover the entire "from"
#           sequence).
chain = import.chain(system.file(package="liftOver", "extdata", "mm9ToMm10.over.chain"))
# chain = import.chain("data/mm9ToMm10.over.chain")

asd <- bulk_bs %>% copy %>% 
  .[, chr := paste0("chr", chr)] %>% .[,end:=start] %>% 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE, ignore.strand = TRUE) %>%
  liftOver(chain) %>%
  as.data.frame() %>%
  setDT() %>% 
  .[, .(seqnames, start, rate, N)] %>%
  setnames(c("seqnames","start"), c("chr","pos")) %>%
  .[, chr := gsub("chr", "", chr)] %>%
  setkey(chr, pos)


fwrite(asd, "/Users/ricard/data/Ficz_2013/BS/rep3_parsed_mm10.txt", quote = F, sep = "\t", row.names = F, col.names = T)
