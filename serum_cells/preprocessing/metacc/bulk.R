met <- fread("zcat < /Users/ricard/data/NMT-seq/met/parsed/met_data.tsv.gz", showProgress=F)
acc <- fread("zcat < /Users/ricard/data/NMT-seq/acc/parsed/acc_data.tsv.gz", showProgress=F)

met_bulk <- met[,.(rate=mean(rate), N=.N),by=c("id","anno")]
acc_bulk <- acc[,.(rate=mean(rate), N=.N),by=c("id","anno")]


fwrite(met_bulk,"/Users/ricard/data/NMT-seq/met/parsed/met_data_pseudobulk.tsv", sep="\t")
fwrite(acc_bulk,"/Users/ricard/data/NMT-seq/acc/parsed/acc_data_pseudobulk.tsv", sep="\t")
