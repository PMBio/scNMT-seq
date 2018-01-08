# Script to calculate binary methylation states from bismark files 

# opts$input_format=1
# chr     pos     rate
# 1       3019021 0
# 1       3027398 100
# 1       3052955 100

# opts$input_format=2
# chr     pos     met_reads	nonmet_reads
# 1       3019021 1	1
# 1       3027398 0	1
# 1       3052955 5	10

library(data.table)
library(purrr)
library(doParallel)

# Define I/O
io <- list()
io$indir <- "/hps/nobackup/stegle/users/ricard/NMT-seq/met/raw/allele-inspecific/unstranded"
io$outdir <- "/hps/nobackup/stegle/users/ricard/NMT-seq/met/raw/allele-inspecific/unstranded/binarised"; dir.create(io$outdir)

# Define options
opts <- list()

opts$input_format <- 1
opts$remove0.5 <- FALSE # if FALSE, 0.5 sites are rounded to 1 (since there are generally more methylated CpG sites than unmethylated)
opts$cores <- 10

# Process samples
filenames <- list.files(io$indir,pattern="(.tsv.gz)$")
registerDoParallel(cores=opts$cores)
invisible(foreach(i=1:length(filenames)) %dopar% {
# for (i in 1:length(filenames)) {
	sample <- sub(".tsv.gz","",filenames[i])
	print(sample)

	# Load data
	data <- fread(paste("gunzip -c", paste(io$indir,filenames[i],sep="/")), sep="\t", verbose=FALSE)

	# Input format 1 (chr,pos,rate)
	if (opts$input_format == 1) {

		# Define column names
		colnames(data) <- c("chr","pos","rate")

		
		# Methylation rate ranges from 0 to 100
		if (max(data$rate) >= 100) {

			# Sanity check
			tmp <- sum((max(data$rate) > 100) | (min(data$rate) < 0))
			if (tmp>0) cat(sprintf("There are %d CpG sites that have methylation rate higher than 100 or lower than 0\n",tmp))

			# Deal with uncertain sites with rate=50
			if (opts$remove0.5) {
				data <- data[rate!=50,]
			} else {
				data[rate==50,rate:=100]
			}

			# Calculate binary methylation status
			data[,status:=round(rate/100)][,rate:=NULL]
		}

		# Methylation rate ranges from 0 to 1
		else if (max(data$rate) >= 1) {

			# Sanity check
			tmp <- sum((max(data$rate) > 1) | (min(data$rate) < 0))
			if (tmp>0) cat(sprintf("There are %d CpG sites that have methylation rate higher than 100 or lower than 0\n",tmp))

			# Deal with uncertain sites with rate=0.5
			if (opts$remove0.5) {
				data <- data[rate!=0.5,]
			} else {
				data[rate==0.5,rate:=1]
			}

			# Calculate binary methylation status
			data[,status:=round(rate)][,rate:=NULL]
		}

	}

	# Input format 2 (chr,pos,met_reads,nonmet_reads)
	if (opts$input_format == 2) {

		# Define column names
		colnames(data) <- c("chr","pos","met_reads","nonmet_reads")

		# Calculate methylation rates
		data[,rate:=(met_reads/nonmet_reads)]

		# Sanity checks
		tmp <- sum((max(data$rate) > 1) | (min(data$rate) < 0))
		if (tmp>0) cat(sprintf("There are %d CpG sites that have methylation rate higher than 100 or lower than 0\n",tmp))

		# Deal with uncertain sites with rate=50
		if (opts$remove0.5) {
			data <- data[rate!=0.5,]
		} else {
			data[rate==0.5,rate:=1]
		}

		# Calculate binary methylation status
		data[,status:=round(rate)][,rate:=NULL]
	}

	# Save results
	fwrite(data, file=sprintf("%s/%s.tsv",io$outdir,sample), sep="\t", showProgress=FALSE, verbose=FALSE, col.names=FALSE)
})
