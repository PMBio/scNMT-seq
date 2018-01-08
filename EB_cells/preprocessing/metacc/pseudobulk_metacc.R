# met <- fread("zcat < /Users/ricard/data/NMT-seq/met/parsed/met_data.tsv.gz", showProgress=F)
# acc <- fread("zcat < /Users/ricard/data/NMT-seq/acc/parsed/acc_data.tsv.gz", showProgress=F)

met <- fread("zcat < /Users/ricard/data/NMT-seq/met/parsed/old/submission1/met_data.tsv.gz", showProgress=F)
acc <- fread("zcat < /Users/ricard/data/NMT-seq/acc/parsed/old/submission1/acc_data.tsv.gz", showProgress=F)

#######
# met[,metreads:=round((rate*weight)/100,0)]
# met_bulk <- met[,.(rate1=mean(rate), rate2=100*sum(metreads)/sum(weight)),by=c("id","anno")] %>% .[,sample:="pseudobulk"]
########

met_bulk <- met[,.(rate=round(mean(rate),2), weight=.N),by=c("id","anno")] %>% .[,sample:="pseudobulk"]
acc_bulk <- acc[,.(rate=round(mean(rate),2), weight=.N),by=c("id","anno")] %>% .[,sample:="pseudobulk"]

# fwrite(met_bulk,"/Users/ricard/data/NMT-seq/met/parsed/met_data_pseudobulk.tsv", sep="\t")
# fwrite(acc_bulk,"/Users/ricard/data/NMT-seq/acc/parsed/acc_data_pseudobulk.tsv", sep="\t")

fwrite(met_bulk,"/Users/ricard/data/NMT-seq/met/parsed/old/submission1/met_data_pseudobulk.tsv", sep="\t")
fwrite(acc_bulk,"/Users/ricard/data/NMT-seq/acc/parsed/old/submission1/acc_data_pseudobulk.tsv", sep="\t")

system("gzip -f /Users/ricard/data/NMT-seq/met/parsed/old/submission1/met_data_pseudobulk.tsv")
system("gzip -f /Users/ricard/data/NMT-seq/acc/parsed/old/submission1/acc_data_pseudobulk.tsv")