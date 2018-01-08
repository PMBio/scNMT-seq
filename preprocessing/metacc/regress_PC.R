
# Script to regress out principal components or known covariates from the methylation data
# In particular, we want to regress out the mean methylation rate

library(data.table)
library(purrr)
# source("/Users/ricard/NOMe-seq/Rutils/stats_utils.R")

# accessibility
# infile <- "/Users/ricard/data/NMT-seq/dna/parsed/allele_unspecific/filt/data.rds"
# outfile <- "/Users/ricard/data/NMT-seq/dna/parsed/allele_unspecific/filt/data_regressed.rds"
# statsfile <- "/Users/ricard/NOMe-seq/stats/samples/GC/sample_stats.txt"

# methylation
infile <- "/Users/ricard/data/NMT-seq/met/parsed/met_data.tsv.gz"
outfile <- "/Users/ricard/data/NMT-seq/met/parsed/met_data_regressed.tsv"
statsfile <- "/Users/ricard/NMT-seq/stats/samples/sample_stats.txt"

# Define annotations to use
# annos <- c("prom_cgi","prom_noncgi","prom_active","prom_inactive","prom_poised","genebody","active_enhancers","super_enhancers","primed_enhancers","CTCF","p300","Nanog","Oct4","IAP")

# Load methylation data
# data <- readRDS(infile) %>% data.table %>% .[,m:=NULL] %>% .[anno %in% annos] %>% .[,anno:=gdata::drop.levels(anno)]
data <- fread(sprintf("zcat < %s",infile))

# Load stats file
stats <- fread(statsfile, header=T)[context=="CG"] %>% setnames("mean","cell_rate") 

#################################
## Regress out known covariate ##
#################################

# How is the regression performed?
# (1) Impute missing values: cell-wise 
# (2) Center the data: since we are computing correlations we don't care about the mean value
# (Q) SHOULD WE USE THE BETA OR M VALUE FOR CORRELATIONS?

# Define the covariate (centered and scaled)
covariate <- (stats$cell_rate - min(stats$cell_rate)) / (max(stats$cell_rate) - min(stats$cell_rate))
names(covariate) <- stats$sample

## Option 1: convert to matrix, impute missing values and do linear model

# Convert from wide to long format
# foo <- data %>% split(.$anno) %>% map(~ dcast(.[,c("sample","id","rate")], sample~id,value.var="rate") %>% as.data.frame %>% tibble::column_to_rownames("sample") %>% as.matrix)
# Center the data
# ....
# Impute missing values
# bar <- foo %>% map(impute,margin=1)
# Do linear model
# baz <- bar %>% map(~ lm(.~covariate) %>% .$residuals)
# Extract residuals

## Option 2: linear model directly to the data.frame

# Filter 
N <- 20
data <- data[,n:=.N,by=c("id","anno")] %>% .[n>=N] %>% .[,n:=NULL]

foo <- data %>% split(.$anno) %>%
  # Center the data
  map(~ .[,mean:=mean(rate),by="id"] %>% .[,rate:=rate-mean] %>% .[,mean:=NULL] %>%
  # Do linear model
  .[,residual:=lm(formula=.SD[,rate]~covariate[.SD[,sample]])[["residuals"]], by=c("id","anno")] %>%
  .[,c("rate","residual"):=.(round(rate,2),round(residual,2))]
  ) %>% rbindlist
  
  

# asd <- foo[,.(l=lm(formula=.SD[,rate]~covariate[.SD[,sample]])[c("coefficients","residuals","fitted.values")]), by="id"]
# asd <- head(foo,n=1000)[,.(l=lm(formula=.SD[,rate]~covariate[.SD[,sample]])[["residuals"]]), by="id"]

# tmp <- foo[id=="ENSMUSG00000000001"]
# tmp1 <- lm(formula=tmp[,rate]~covariate[tmp[,sample]])

# Save data
# saveRDS(foo,outfile)
fwrite(foo, outfile, sep="\t")

####################
## Regress out PC ##
####################
# 
# foo <- tmp %>% split(.$anno) %>% map(~ dcast(.[,c("sample","id","rate")], sample~id,value.var="rate") %>% as.data.frame %>% tibble::column_to_rownames("sample") %>% as.matrix)
# bar <- foo %>% map(impute,margin=1)
# baz <- bar %>% map(prcomp,scale=F,center=T)
# 
# p <- baz[[1]]
# plot_pca_val(p$sdev**2/sum(p$sdev**2))
# abs(cor(stats$cell_rate,p$x))
# 
# stop()

# Regress out PC1
# tmp <- t()
# plot(as.vector(tmp),as.vector(p$rotation[,1]%*%t(matrix(p$x[,1]))))
# tmp <- tmp - t(p$rotation[,1] %*% t(matrix(p$x[,1])))

# Check by recomputing PCA
# p <- prcomp(tmp,scale=F,center=T)
# plot_pca_val(p$sdev**2/sum(p$sdev**2))
# abs(cor(stats$cell_rate,p$x))


##########
## Test ##
##########

# test_id <- "active_enhancers_3414"
# to.plot <- tmp[id==test_id] %>% drop.levels()
# asd <- to.plot %>% setorder(cell_rate) %>% .$sample
# to.plot1 <- melt(to.plot, measure.vars=c("rate","cell_rate"), variable.name="rate_type", value.name="value")%>% .[,sample:=factor(sample,levels=asd)]
#  
# ggplot(to.plot1, aes(x=sample, y=value, fill=rate_type)) +
#   geom_bar(stat="identity", position="dodge")
#   
#   

