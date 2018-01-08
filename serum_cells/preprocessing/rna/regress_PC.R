
#############################################################################################
## Script to regress out principal components or known covariates from the expression data ##
#############################################################################################


library(data.table)
library(purrr)
library(sceset)
source("/Users/ricard/NOMe-seq/Rutils/stats_utils.R")

####################
## Define options ##
####################

# I/O
io <- list()
io$infile <- "/Users/ricard/data/NMT-seq/rna/parsed/allele_unspecific/sceset.rds"
io$outfile <- "/Users/ricard/data/NMT-seq/rna/parsed/allele_unspecific/sceset_regressed.rds"

# Options
opts <- list()
opts$center <- T # Center all genes to zero expression

###############
## Load data ##
###############

# Load expression data
sceset <- readRDS(io$infile)

###############
## Parse data ##
###############

# Keep N sites with largest variability
# N <- "all"
# if (N=="all") {
#   genes.keep <- apply(exprs(sceset),1,var) %>% sort %>% names
# } else {
#   genes.keep <- apply(exprs(sceset),1,var) %>% sort %>% tail(n=N) %>% names
# }
# sceset_filt <- sceset[genes.keep,]

# Center the genes
if (opts$center == T) {
  exprs(sceset) <- sweep(x=exprs(sceset), MARGIN=1, STATS=apply(exprs(sceset),1,mean), FUN="-")
}
  
#####################
## Regress out PCs ##
#####################

# Do PCA
p <- prcomp(t(exprs(sceset)), scale=F, center=T)
plot_pca_val(p$sdev**2/sum(p$sdev**2))

# Check correlation of PCs to known drivers of variation
abs(cor(sceset$pct_dropout,p$x))
abs(cor(sceset$total_features,p$x))

# Regress out PC1
tmp <- t(exprs(sceset))
# plot(as.vector(tmp),as.vector(p$rotation[,1]%*%t(matrix(p$x[,1]))))
tmp <- tmp - t(p$rotation[,1] %*% t(matrix(p$x[,1])))

# Check by recomputing PCA
p <- prcomp(tmp,scale=F,center=T)
plot_pca_val(p$sdev**2/sum(p$sdev**2))
abs(cor(sceset$pct_dropout,p$x))
abs(cor(sceset$total_features,p$x))

############################
## Regress out covariates ##
############################

# Fit the linear model
covariate <- (sceset$pct_dropout - min(sceset$pct_dropout)) / (max(sceset$pct_dropout) - min(sceset$pct_dropout))
r <- lm(t(exprs(sceset)) ~ covariate)
# View(t(r$coefficients))

# # Plot to test
# gene="Arfgef1"
# r$coefficients[,gene]
# predicted <- as.matrix(data.frame(1,covariate)) %*% as.matrix(r$coefficients[,gene])
# true <- exprs(sceset)[gene,]
# plot(covariate,true)
# # plot(predicted,true)
# abline(coef=r$coefficients[,gene])

# Regress out 
# regressed <- exprs(sceset) - t(as.matrix(data.frame(1,covariate)) %*% r$coefficients)
exprs(sceset) <- t(r$residuals)

# Check again if the covariate explains variation
# r <- lm(t(exprs(sceset)) ~ covariate)

###############
## Save data ##
###############

sceset <- calculateQCMetrics(sceset)

saveRDS(sceset,file=out)
