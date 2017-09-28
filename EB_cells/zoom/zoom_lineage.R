
##########################
## Zoom plot by lineage ##
##########################

suppressMessages(library(scater))
suppressMessages(library(data.table))
suppressMessages(library(purrr))
suppressMessages(library(ggplot2))
suppressMessages(library(argparse))

source("/Users/ricard/NMT-seq/zoom/andreas/utils_region.R")

####################
## Define options ##
####################

## Options ##
opts <- list()
opts$gene <- "Prtg"
opts$window <- 11000
opts$slide <- 1000
opts$up <- 20000
opts$down <- 20000


## I/O ##
io <- list()
io$basedir <- "/Users/ricard/data/NMT-seq_EB"
io$plot.outdir <- "/Users/ricard/NMT-seq/rebuttal/EB/zoom/out"
io$genes.infile <- "/Users/ricard/data/ensembl/mouse/v87/BioMart/mRNA/Mmusculus_genes_BioMart.87.txt"
io$in.sample_metadata <- paste(io$basedir,"sample_sheet.csv",sep="/")
io$met_dir <- paste0(io$basedir,"/met/raw/merged/binarised/")
io$acc_dir <- paste0(io$basedir,"/acc/raw/merged/binarised/")
io$rna_file <- paste0(io$basedir,"/rna/parsed/sceset.rds")
io$features_dir <- paste0(io$basedir,"/features/filt")
# io$plot.outfile <- sprintf("%s/%s.pdf", io$plot.outdir, opts$gene)

## Options ##
opts$features <- c("CGI","prom_2000_0","exons","active_enhancers","super_enhancers","Nanog","Oct4","p300","CTCF")
opts$cells <- fread(io$in.sample_metadata, header=T) %>% .[pass_metQC==T & pass_rnaQC==T & pass_accQC==T,sample]

###############
## Load data ##
###############

## Load gene information
gene <- fread(sprintf(io$genes.infile))[,c(1,2,3,4,7)] %>%
  .[,chr:=sub("chr","",chr)] %>%
  .[symbol==opts$gene]

opts$chr <- gene$chr
opts$start <- gene$start-opts$up
opts$end <- gene$end+opts$down

# Methylation and Accessibility 
filename <- sprintf(paste0(io$plot.outdir, "/%s_chr%s_%s-%s.rds"),opts$gene,opts$chr,opts$start,opts$end)
if (file.exists(filename)) {
  data <- readRDS(filename)
} else {
  data <- read_genomic_windows(io,opts)
  saveRDS(data, file=filename)
}
# data$met[,rate:=rate*100]

# Expression
rna <- exprs(readRDS(io$rna_file))[opts$gene,]
data$rna <- data.table(sample=names(rna), expr=rna)

# Genomic contexts
files <- paste(io$features_dir,paste0(opts$features,".bed"),sep="/")
anno_list <- files %>% map(~ fread(.x, sep="\t", verbose=F, stringsAsFactors=F)[,c(1,2,3,4,5,6)])
anno_df <- rbindlist(anno_list) %>% setnames(c("chr","start","end","strand","ens_id","anno")) %>% 
  .[,anno:=as.factor(anno)] %>% .[,chr:=as.factor(sub("chr","",chr))] %>% setkey("chr","start","end") %>%
  foverlaps(data.table(chr=opts$chr,start=opts$start,end=opts$end) %>% setkey("chr","start","end"), nomatch=0) %>%
  .[,c("chr","i.start","i.end","anno")] %>% setnames(c("i.start","i.end"),c("start","end"))
anno_df$start <- ifelse(anno_df$start<gene$start-opts$up, gene$start-opts$up, anno_df$start)
anno_df$end <- ifelse(anno_df$end>gene$end+opts$down, gene$end+opts$down, anno_df$end)

##########################
## Perform calculations ##
##########################

# Calculate running means, variances and correlations
tmp <- seq(from=opts$start, to=opts$end-opts$window, by=opts$slide)
foo <- data.table(
  window_center=tmp+(opts$window/2),
  chr=opts$chr,
  start=tmp,
  end=tmp+opts$window
)
bar <- rbind(data$met %>% copy %>% .[,context:="CG"], data$acc %>% copy %>% .[,context:="GC"]) %>%
  .[,c("start","end","pos") := list(pos,pos,NULL)]
to.plot <- foverlaps( foo%>%setkey("chr","start","end"), bar %>%setkey("chr","start","end")) %>%
  .[,c("chr","i.start","i.end","start","end"):=NULL] %>%
  .[,.(rate=mean(rate), N=.N),by=c("sample","window_center","context")]


##############################
## Merge data with metadata ##
##############################

sample_metadata <- fread(io$in.sample_metadata, header=T) %>% .[sample%in%opts$cells]
to.plot <- merge(to.plot,sample_metadata[,c("sample","lineage")], by="sample")

####################
## Generate plots ##
####################

p <- list()

## Plot profile ##

asd <- to.plot[,.(mean=mean(rate)), by=c("window_center","context")] %>%
  .[,.(min=min(mean), max=max(mean)),by="context"] %>% setkey(context)

# f <- function(x) { return(data.frame(y=mean(x), ymin=mean(x)-sd(x)/2, ymax=mean(x)+sd(x)/2)) }
mean_f <- function(x,context) { return( 100 * (mean(x) - asd[context,min]) / ( asd[context,max] - asd[context,min])) }
all_f <- function(x,context) { return(
  data.frame(
    y = 100 * (mean(x) - asd[context,min]) / ( asd[context,max] - asd[context,min]),
    ymin = 100 * (mean(x) - asd[context,min]) / ( asd[context,max] - asd[context,min]) - sd(x)/2.0,
    ymax = 100 * (mean(x) - asd[context,min]) / ( asd[context,max] - asd[context,min]) + sd(x)/2.0
  )
) }

# p$profile1 <- ggplot(to.plot[context=="GC" & lineage=="Pluripotent"], aes(x=window_center, y=rate, color=context, fill=context)) +
#   stat_summary(fun.y=mean_f, fun.args="GC", geom="line", size=1.0, data=to.plot[context=="GC" & lineage=="Pluripotent"]) +
#   stat_summary(fun.data=all_f, fun.args="GC", geom="smooth", alpha=0.3, linetype=0, fill="#00BFC4", data=to.plot[context=="GC" & lineage=="Pluripotent"]) +
#   stat_summary(fun.y=mean_f, fun.args="CG", geom="line", size=1.0, data=to.plot[context=="CG" & lineage=="Pluripotent"]) +
#   stat_summary(fun.data=all_f, fun.args="CG", geom="smooth", alpha=0.3, linetype=0, fill="#F8766D", data=to.plot[context=="CG" & lineage=="Pluripotent"]) +
#   scale_y_continuous(limits=c(-20,120), breaks=c(0,20,40,60,80,100)) +
#   scale_x_continuous(limits=c(opts$start,opts$end)) +
#   scale_colour_manual(labels=c("CG methylation", "GC accessibility"), values=c("#F8766D","#00BFC4")) +
#   guides(fill=FALSE) +
#   xlab("") + ylab("methylation/accessibility scaled rate") +
#   theme(
#     axis.text.x = element_blank(),
#     axis.text.y = element_blank(),
#     axis.title.y = element_blank(),
#     axis.line = element_blank(),
#     axis.ticks.x = element_blank(),
#     axis.ticks.y = element_line(size=rel(1.0), color="black"),
#     legend.position = "top",
#     legend.direction = "horizontal",
#     legend.key = element_blank(),
#     legend.margin = margin(t=0, r=0, b=-25, l=0, unit="pt"),
#     legend.title = element_blank(),
#     legend.text = element_text(size=13),
#     panel.border=element_blank(),
#     panel.grid.major=element_blank(),
#     panel.grid.minor=element_blank(),
#     panel.background=element_blank()
#   )
# p$profile1

# p$profile2 <- ggplot(to.plot[context=="GC" & lineage=="Differentiated"], aes(x=window_center, y=rate, color=context, fill=context)) +
#   stat_summary(fun.y=mean_f, fun.args="GC", geom="line", size=1.0, data=to.plot[context=="GC" & lineage=="Differentiated"]) +
#   stat_summary(fun.data=all_f, fun.args="GC", geom="smooth", alpha=0.3, linetype=0, fill="#00BFC4", data=to.plot[context=="GC" & lineage=="Differentiated"]) +
#   stat_summary(fun.y=mean_f, fun.args="CG", geom="line", size=1.0, data=to.plot[context=="CG" & lineage=="Differentiated"]) +
#   stat_summary(fun.data=all_f, fun.args="CG", geom="smooth", alpha=0.3, linetype=0, fill="#F8766D", data=to.plot[context=="CG" & lineage=="Differentiated"]) +
#   scale_y_continuous(limits=c(-20,120), breaks=c(0,20,40,60,80,100)) +
#   scale_x_continuous(limits=c(opts$start,opts$end)) +
#   scale_colour_manual(labels=c("CG methylation", "GC accessibility"), values=c("#F8766D","#00BFC4")) +
#   guides(fill=FALSE) +
#   xlab("") + ylab("methylation/accessibility scaled rate") +
#   theme(
#     axis.text.x = element_blank(),
#     axis.text.y = element_blank(),
#     axis.title.y = element_blank(),
#     axis.line = element_blank(),
#     axis.ticks.x = element_blank(),
#     axis.ticks.y = element_line(size=rel(1.0), color="black"),
#     legend.position = "top",
#     legend.direction = "horizontal",
#     legend.key = element_blank(),
#     legend.margin = margin(t=0, r=0, b=-25, l=0, unit="pt"),
#     legend.title = element_blank(),
#     legend.text = element_text(size=13),
#     panel.border=element_blank(),
#     panel.grid.major=element_blank(),
#     panel.grid.minor=element_blank(),
#     panel.background=element_blank()
#   )
# p$profile2


p$acc_profile <- ggplot(to.plot, aes(x=window_center, y=rate, color=lineage, fill=lineage)) +
  stat_summary(fun.y=mean, geom="line", size=0.75, data=to.plot[context=="GC" & lineage=="Pluripotent"]) +
  stat_summary(fun.y=mean, geom="line", size=0.75, data=to.plot[context=="GC" & lineage=="Differentiated"]) +
  stat_summary(fun.data=mean_se, geom="smooth", alpha=0.2, linetype=0, fill="dodgerblue1", data=to.plot[context=="GC" & lineage=="Pluripotent"]) +
  stat_summary(fun.data=mean_se, geom="smooth", alpha=0.2, linetype=0, fill="dodgerblue4", data=to.plot[context=="GC" & lineage=="Differentiated"]) +
  # scale_y_continuous(limits=c(-20,120), breaks=c(0,20,40,60,80,100)) +
  # scale_x_continuous(limits=c(opts$start,opts$end)) +
  scale_colour_manual(labels=c("Differentiated cells","Pluripotent cells"), values=c("dodgerblue4","dodgerblue1")) +
  xlab("") + ylab("Accessibility rate") +
  guides(fill=FALSE) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=rel(1.3), color="black"),
    axis.title.y = element_text(size=rel(1.5), color="black"),
    # axis.line = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(size=rel(1.0), color="black"),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.key = element_blank(),
    legend.margin = margin(t=0, r=0, b=-15, l=0, unit="pt"),
    legend.title = element_blank(),
    legend.text = element_text(size=13),
    panel.border=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank()
  )
# p$acc_profile

p$met_profile <- ggplot(to.plot, aes(x=window_center, y=rate, color=lineage, fill=lineage)) +
  stat_summary(fun.y=mean, geom="line", size=0.75, data=to.plot[context=="CG" & lineage=="Pluripotent"]) +
  stat_summary(fun.y=mean, geom="line", size=0.75, data=to.plot[context=="CG" & lineage=="Differentiated"]) +
  stat_summary(fun.data=mean_se, geom="smooth", alpha=0.2, linetype=0, fill="orangered1", data=to.plot[context=="CG" & lineage=="Pluripotent"]) +
  stat_summary(fun.data=mean_se, geom="smooth", alpha=0.2, linetype=0, fill="orangered4", data=to.plot[context=="CG" & lineage=="Differentiated"]) +
  scale_y_continuous(limits=c(-10,110), breaks=c(0,25,50,75,100)) +
  # scale_x_continuous(limits=c(opts$start,opts$end)) +
  scale_colour_manual(labels=c("Differentiated cells", "Pluripotent cells"), values=c("orangered4","orangered1")) +
  xlab("") + ylab("Methylation rate") +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=rel(1.3), color="black"),
    axis.title.y = element_text(size=rel(1.5), color="black"),
    # axis.line = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(size=rel(1.0), color="black"),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.key = element_blank(),
    legend.margin = margin(t=0, r=0, b=-15, l=0, unit="pt"),
    legend.title = element_blank(),
    legend.text = element_text(size=13),
    panel.border=element_blank(),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_blank()
  )
# p$met_profile

# 
# to.plot_relative <- to.plot[,.(rate=mean(.SD[lineage=="Pluripotent",rate]) - mean(.SD[lineage=="Differentiated",rate])), by=c("window_center","context")]
# p$relative_profile <- ggplot(to.plot_relative, aes(x=window_center, y=rate, color=context, fill=context)) +
#   geom_line(aes(color=context)) +
#   scale_colour_manual(labels=c("CG methylation", "GC accessibility"), values=c("#F8766D","#00BFC4")) +
#   guides(fill=FALSE) +
#   xlab("") + ylab("methylation/accessibility scaled rate") +
#   theme(
#     legend.position = "top",
#     legend.direction = "horizontal",
#     legend.key = element_blank(),
#     legend.margin = margin(t=0, r=0, b=-25, l=0, unit="pt"),
#     legend.title = element_blank(),
#     legend.text = element_text(size=13),
#     panel.border=element_blank(),
#     panel.grid.major=element_blank(),
#     panel.grid.minor=element_blank(),
#     panel.background=element_blank()
#   )
# p$relative_profile


## Plot correlation ##
cor <- to.plot %>% 
  dcast(sample+window_center~context, value.var="rate") %>%
  merge(data$rna,by="sample") %>%
  .[,.(r_metrna=cor(x=CG, y=expr, method=opts$cor_method, use="pairwise.complete.obs"),
       r_metacc=cor(x=CG, y=GC, method=opts$cor_method, use="pairwise.complete.obs"),
       r_accrna=cor(x=GC, y=expr, method=opts$cor_method, use="pairwise.complete.obs")), by=window_center]  %>%
  melt(id.vars="window_center",variable.name="r_type")

p$cor <- ggplot(cor,aes(x=window_center,y=value)) +
  geom_line(aes(color=r_type), size=0.6) +
  scale_x_continuous(limits=c(opts$start,opts$end)) +
  scale_y_continuous(limits=c(-1.0,1.0), breaks=c(-0.5,0,0.5)) +
  # geom_hline(yintercept=0, linetype="dashed", alpha=0.75, size=0.75) +
  geom_segment(x=opts$start, xend=opts$end, y=0, yend=0, color="black", size=0.6,  linetype="dashed") +
  ylab("Correlation") +
  scale_colour_discrete(labels=c("Met/Exp","Met/Acc","Acc/Expr")) +
  # scale_colour_discrete(labels=c("Methylation/Expression","Methylation/Accessibility","Accessibility/Expression")) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=rel(1.3), color="black"),
    axis.title.y = element_text(size=rel(1.5), color="black"),
    axis.title.x = element_blank(),
    # axis.title.y = element_text(size=13, margin=margin(0,10,0,0)),
    axis.line = element_blank(),
    # axis.line = element_line(size=rel(1.0)),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(size=rel(1.1), color="black"),
    legend.key = element_blank(),
    legend.position = "top",
    # legend.position = c(0.5,0.1),
    legend.direction = "horizontal",
    # legend.key.width = unit(2.0,"line"),
    # legend.key.height = unit(2.0,"line"),
    legend.margin = margin(t=0, r=0, b=-10, l=0, unit="pt"),
    legend.title = element_blank(),
    legend.text = element_text(size=13),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()
  )
# p$cor


## Plot annotation track ##
# p$anno <- ggplot(anno_df) +
#   geom_rect(aes(xmin=start, xmax=end, ymin=as.numeric(anno)-0.20, ymax=as.numeric(anno)+0.20, fill=anno), alpha=1.0) +
#   scale_x_continuous(limits=c(opts$start,opts$end)) +
#   theme(
#     # plot.margin = unit(c(t=1,r=1,b=1,l=1), "cm"),
#     plot.title = element_blank(),
#     axis.text.y=element_blank(),
#     axis.text.x=element_text(color="black", size=rel(1.0)),
#     axis.ticks.x = element_line(color="black", size=rel(1.1)),
#     axis.ticks.y = element_blank(),
#     axis.title=element_blank(),
#     axis.line = element_line(color="black", size=rel(0.7)),
#     legend.key = element_blank(),
#     legend.position = "bottom",
#     legend.direction = "horizontal",
#     # legend.key.size= unit(0.5, "cm"),
#     legend.key.width=unit(1.0,"line"),
#     legend.key.height=unit(1.0,"line"),
#     # legend.margin = margin(t=10, r=0, b=0, l=0, unit="pt"),
#     legend.title = element_blank(),
#     legend.text = element_text(size=7),
#     panel.border=element_blank(),
#     panel.grid.major=element_blank(),
#     panel.grid.minor=element_blank(),
#     panel.background=element_blank()
#   )

## Join all plots
# final_plot <- cowplot::plot_grid(p$cor, p$acc_profile, p$met_profile, p$anno, align="v", nrow=4, rel_heights=c(1/6,1/5,1/5,1/6))
final_plot <- cowplot::plot_grid(p$cor, p$acc_profile, p$met_profile, align="v", nrow=3, rel_heights=c(1/6,1/5,1/5))

###############
## Save plot ##
###############

final_plot

io$plot.outfile <- paste0(io$plot.outdir,"/",opts$gene,".pdf")
pdf(file=io$plot.outfile, width = 5, height=7)
print(final_plot)
dev.off()


