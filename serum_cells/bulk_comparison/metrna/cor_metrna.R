###############################################################
## Correlation between bulk methylation and RNA expression ##
###############################################################

############################
## Define I/O and options ##
############################

library(data.table)
library(purrr)
library(biomaRt)

io <- list()
io$rna <- "/Users/ricard/data/encode/RNA/Encode_E14_totalRNAseq_gene_quantifications.tsv.gz"
io$met <- "/Users/ricard/data/Ficz_2013/BS/parsed/met_data.tsv.gz"

###############
## Load data ##
###############

rna = fread(paste0("zcat < ",io$rna))
met = fread(paste0("zcat < ",io$met))

mart = useMart("ENSEMBL_MART_ENSEMBL") %>% useDataset("mmusculus_gene_ensembl", .)
genes_mm10 = getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "external_gene_name"), 
                   filters=c("biotype"), values=list("protein_coding"), mart = mart) %>% setDT() %>% setnames(c("gene_id", "chr", "start", "end", "gene"))

################
## Parse data ##
################

met <- met  %>% .[,m:=log2(((rate/100)+0.01)/(1-(rate/100)+0.01))]

# gene ids seem to be separated by transcript?
rna = tidyr::separate(rna, gene_id, c("gene_id", "transcript"), sep = "\\.", fill = "right")


rna <- rna[, .(gene_id, FPKM)] %>%
  setkey(gene_id) %>%
  merge(genes_mm10 %>% setkey(gene_id)) %>%
  .[, log_FPKM := log2(FPKM + 1)] %>%
  .[, start := start - 10000] %>%
  .[, end := end + 10000]


######################################
## Merge RNA and methylation data ##
######################################

metrna = setkey(met, chr, start, end) %>%
  foverlaps(rna %>% setkey(chr, start, end), nomatch = 0)


##########################################
## Correlate RNA and methylation data ##
##########################################

cor = metrna[, .(r = round(cor(log_FPKM, m, method="pearson"),3)), anno]

##################
## Save results ##
##################

fwrite(cor, "/Users/ricard/NMT-seq/bulk_comparison/metrna/out/cor_metrna.txt", quote=F, col.names = T, row.names = F, sep="\t")

