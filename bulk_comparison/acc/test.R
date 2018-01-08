############################
## Define I/O and options ##
############################

library(data.table)
library(purrr)

io <- list()
io$acc <- "/Users/ricard/data/encode/DNAse-seq/parsed/encode_DNAseq.tsv.gz"

###############
## Load data ##
###############

acc = fread(paste0("zcat < ",io$acc))


###############
## Do stuff ##
###############

acc[,.(rate=m)]
p <- ggplot(accrna[anno=="prom_50_50"], aes(x=dnase_log,y=log_FPKM)) +
  geom_point(size=0.4, alpha=0.5) +
  stat_smooth(method="lm", color="blue", alpha=0.5) +
  # facet_wrap(~sample) +
  # ggtitle("Active Enhancers (ENCODE)") +
  xlab("DNAse-seq log normalised counts") + ylab("RNA-seq log normalised counts") +
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
p
