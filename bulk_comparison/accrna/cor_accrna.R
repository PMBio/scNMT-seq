###############################################################
## Correlation between bulk accessibility and RNA expression ##
###############################################################

############################
## Define I/O and options ##
############################

library(data.table)
library(purrr)
library(biomaRt)

io <- list()
io$rna <- "/Users/ricard/data/encode/RNA/Encode_E14_totalRNAseq_gene_quantifications.tsv.gz"
io$acc <- "/Users/ricard/data/encode/DNAse-seq/parsed/encode_DNAseq.tsv.gz"

###############
## Load data ##
###############

rna = fread(paste0("zcat < ",io$rna))

acc = fread(paste0("zcat < ",io$acc))

mart = useMart("ENSEMBL_MART_ENSEMBL") %>% useDataset("mmusculus_gene_ensembl", .)
genes_mm10 = getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "external_gene_name"), 
                   filters=c("biotype"), values=list("protein_coding"), mart = mart) %>% setDT() %>% setnames(c("gene_id", "chr", "start", "end", "gene"))

################
## Parse data ##
################

# gene ids seem to be separated by transcript?
rna = tidyr::separate(rna, gene_id, c("gene_id", "transcript"), sep = "\\.", fill = "right")

# Calculate log FPKM
rna <- rna[, .(gene_id, FPKM)] %>%
  setkey(gene_id) %>%
  merge(genes_mm10 %>% setkey(gene_id)) %>%
  .[, log_FPKM := log2(FPKM + 1)] %>%
  .[, start := start - 10000] %>%
  .[, end := end + 10000]


######################################
## Merge RNA and accessibility data ##
######################################

accrna = setkey(acc, chr, start, end) %>%
  foverlaps(rna %>% setkey(chr, start, end), nomatch = 0)


##########################################
## Correlate RNA and accessibility data ##
##########################################

cor = accrna[, .(r = round(cor(log_FPKM, dnase_log, method="pearson"),3)), anno]

##################
## Scatterplots ##
##################


# p <- ggplot(accrna[anno=="prom_50_50" & dnase_log>-5 & dnase_log<2.5], aes(x=dnase_log,y=log_FPKM)) +
# p <- ggplot(accrna[anno=="activeEnhancers_ENCODE" & dnase_log>-7.5 & dnase_log<2.5], aes(x=dnase_log,y=log_FPKM)) +
p <- ggplot(accrna[anno=="activeEnhancers_Creyghton" & dnase_log>-7.5 & dnase_log<(-1)], aes(x=dnase_log,y=log_FPKM)) +
  # geom_point(size=0.4, alpha=0.5) +
  stat_binhex(bins=250) +
  # geom_bin2d(bins=250) +
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

pdf("/Users/ricard/NMT-seq/bulk_comparison/accrna/out/bulk_enhancerCreyghton.pdf")
print(p)
dev.off()

##################
## Save results ##
##################

fwrite(cor, "/Users/ricard/NMT-seq/bulk_comparison/accrna/out/cor_accrna.txt", quote=F, col.names = T, row.names = F, sep="\t")
