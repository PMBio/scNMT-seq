suppressMessages(library(data.table))
suppressMessages(library(purrr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggbio))
source("/Users/ricard/NMT-seq/zoom/andreas/utils_region.R")

####################
## Define options ##
####################


## I/O ##
io <- list()
io$basedir <- "/Users/ricard/data/NMT-seq_EB"
io$genes.infile <- "/Users/ricard/data/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt"
io$features_dir <- paste0(io$basedir,"/features/filt")

## Options ##

opts <- list()
opts$gene <- "Prtg"
opts$up <- 20000
opts$down <- 20000
# opts$features <- c("prom_2000_0","exons","active_enhancers","Nanog","p300")
opts$features <- c("prom_2000_0","exons","super_enhancers","p300")


###############
## Load data ##
###############

# Load gene information
gene <- fread(sprintf(io$genes.infile))[,c(1,2,3,4,7)] %>%
  .[,chr:=sub("chr","",chr)] %>%
  .[symbol==opts$gene]

opts$chr <- gene$chr
opts$start <- gene$start-opts$up
opts$end <- gene$end+opts$down

# Load Genomic contexts
files <- paste(io$features_dir,paste0(opts$features,".bed"),sep="/")
anno_list <- files %>% map(~ fread(.x, sep="\t", verbose=F, stringsAsFactors=F)[,c(1,2,3,4,5,6)])
anno_df <- rbindlist(anno_list) %>% setnames(c("chr","start","end","strand","ens_id","anno")) %>% 
  .[,anno:=as.factor(anno)] %>% .[,chr:=as.factor(sub("chr","",chr))] %>% setkey("chr","start","end") %>%
  foverlaps(data.table(chr=opts$chr,start=opts$start,end=opts$end) %>% setkey("chr","start","end"), nomatch=0) %>%
  .[,c("chr","i.start","i.end","anno")] %>% setnames(c("i.start","i.end"),c("start","end"))
anno_df$start <- ifelse(anno_df$start<gene$start-opts$up, gene$start-opts$up, anno_df$start)
anno_df$end <- ifelse(anno_df$end>gene$end+opts$down, gene$end+opts$down, anno_df$end)

##########
## Plot ##
##########



# opts$features <- c("CGI","prom_2000_0","exons","active_enhancers","Nanog","Oct4","p300")


asd <- GenomicRanges::makeGRangesFromDataFrame(anno_df, keep.extra.columns=T)
p1 <- autoplot(asd[asd$anno=="exons"], fill="black", geom="alignment")
p2 <- autoplot(asd[asd$anno=="prom_2000_0"], fill="red")
p3 <- autoplot(asd[asd$anno=="super_enhancers"], fill="green")
p4 <- autoplot(asd[asd$anno=="p300"], fill="purple")

# p$anno <- p1 + theme_bw() + scale_x_continuous(limits=c(opts$window_start,opts$window_end))
#
p <- tracks('Exons'=p1, 'Promoter'= p2, 'Super enhancers'=p3, 'p300'=p4, 
            xlab="", main="", 
            scale.height=unit(5.0,"lines"),
            label.bg.fill="white", label.text.color="white", label.text.angle=0, label.text.cex=0.0) +
  # scale_x_continuous(limits=c(opts$window_start,opts$window_end)) +
  theme_bw() +
  theme(
    # plot.margin = unit(c(t=1,r=1,b=1,l=5), "cm"),
    axis.text = element_blank(),
    axis.title=element_blank(),
    # axis.line = element_blank(),
    axis.ticks = element_blank(),
    # panel.border=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank()
    # panel.background=element_blank()    
  )
print(p)

pdf("/Users/ricard/NMT-seq/rebuttal/EB/zoom/out/Prtg_track.pdf", height=3, width=5)
print(p)
dev.off()
