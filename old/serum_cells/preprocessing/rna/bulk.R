# Create a bulk genome

library(scater)

io <- list()
io$infile <- "/Users/ricard/data/NMT-seq/rna/parsed/allele_unspecific/sceset.rds"
io$outfile <- "/Users/ricard/data/NMT-seq/rna/parsed/allele_unspecific/sceset_bulk.rds"

sce <- readRDS(io$infile)

sce_bulk <- newSCESet(exprsData = as.matrix(apply(exprs(sce),1,mean),nc=1), 
                      # countData = as.matrix(apply(counts(sce),1,sum),nc=1), 
                      phenoData = NULL,
                      featureData = featureData(sce))

saveRDS(sce_bulk,io$outfile)
