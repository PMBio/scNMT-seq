#############################################################
## Correlation between bulk methylation and RNA expression ##
#############################################################

############################
## Define I/O and options ##
############################

io <- list()
io$met <- "/Users/ricard/data/Ficz_2013/BS/parsed/met_data.tsv.gz"
io$acc <- "/Users/ricard/data/encode/DNAse-seq/parsed/encode_DNAseq.tsv.gz"
  

#########################
## Load and parse data ##
#########################

acc = fread(paste0("zcat < ",io$acc))  %>% .[,c("dnase","dnase_norm","dnase_log","id","anno")]
met = fread(paste0("zcat < ",io$met)) %>% .[,.(rate=mean(rate)), by=c("id","anno")] %>%
  .[,m:=log2(((rate/100)+0.01)/(1-(rate/100)+0.01))]


#########################################
## Merge Methylation and accessibility ##
#########################################

metacc <- merge(met,acc,by=c("id","anno"))# %>% .[dnase_norm<10]


#############################################
## Correlate Methylation and Accessibility ##
#############################################

cor = metacc[, .(r = round(cor(rate, dnase_log, method="pearson"),3)), anno]


##################
## Save results ##
##################

fwrite(cor, "/Users/ricard/NMT-seq/bulk_comparison/metacc/out/cor_metacc.txt", quote=F, col.names = T, row.names = F, sep="\t")

