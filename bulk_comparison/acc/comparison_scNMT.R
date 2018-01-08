############################
## Define I/O and options ##
############################

library(data.table)
library(purrr)

io <- list()
io$acc_dnaseq <- "/Users/ricard/data/encode/DNAse-seq/parsed/encode_DNAseq.tsv.gz"
io$acc_scnmt <- "/Users/ricard/data/NMT-seq_EB/acc/parsed/acc_data.tsv.gz"

###############
## Load data ##
###############

acc_dnaseq = fread(paste0("zcat < ",io$acc_dnaseq))
acc_scnmt = fread(paste0("zcat < ",io$acc_scnmt))

###########
## Parse ##
###########

# Pseudobulk scNMT
acc_scnmt_pseudobulk <- acc_scnmt[,.(rate=round(mean(rate),2)), by=c("id","anno")] 

###########
## Merge ##
###########

acc <- merge(
  acc_dnaseq,
  acc_scnmt_pseudobulk,
  by=c("id","anno")
)

acc[,m:=log2(((rate/100)+0.01)/(1-(rate/100)+0.01))]

######################
## Plot correlation ##
######################

annos <- c("DHS","p300","prom_2000_2000")
for (n in annos) {
  print(paste0(n,": ",cor(acc[anno==n,rate],acc[anno==n,dnase_log],method="spearman")))
  p <- ggplot(acc[anno==n], aes(x=dnase_log,y=rate)) +
    geom_point(size=0.4, alpha=0.5) +
    stat_smooth(method="lm", color="blue", alpha=0.5) +
    ggtitle(n) +
    xlab("DNAse-seq log normalised counts") + ylab("Rate") +
    theme(
      plot.title = element_text(size=17, hjust=0.5, margin=margin(0,0,20,0)),
      axis.title.y = element_text(colour="black", size=15, vjust=1.5),
      axis.title.x = element_text(colour="black", size=15, vjust=1.5, margin=margin(10,0,0,0)),
      axis.text.x = element_text(colour="black",size=rel(1.0)),
      axis.text.y = element_text(colour="black",size=rel(1.0)),
      # axis.line = element_line(colour="black", size=rel(0.7)),
      # axis.ticks.x = element_line(colour="black", size=rel(0.8)),
      # axis.ticks.y = element_blank(),
      legend.position="none",
      panel.background = element_blank(),
      panel.border = element_blank(),
      panel.grid = element_blank()
    )
  print(p)
}
