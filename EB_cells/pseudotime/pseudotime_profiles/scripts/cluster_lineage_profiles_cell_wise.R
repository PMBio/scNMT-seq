# ------------------------------------------
# Set working directory and load libraries
# ------------------------------------------
if (interactive()){
  cur.dir <- dirname(parent.frame(2)$ofile)
  setwd(cur.dir)
}

suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(BPRMeth))
suppressPackageStartupMessages(library(pracma))
suppressPackageStartupMessages(library(data.table))

# Initialize parameters
io            <- list(type = "acc")
io$base_dir   <- "/home/cakapourani/Documents/Code/datasets/nmt-seq-data-eb"
# io$base_dir <- "/Users/ricard/data/NMT-seq"
# io$base_dir <- "/Users/kapouranis/Documents/Andreas/Code/datasets/nmt-seq-data"
io$data_dir   <- paste0(io$base_dir, "/", io$type, "/parsed/profiles/")
# Load saved objects for accessibility data
load(paste0(io$data_dir, io$type, "_400bp.RData"))
io$data_dir   <- paste0(io$base_dir, "/", io$type, "/parsed/profiles/")

annos              <- "prom_200_200"
opts$basis_prof    <- create_rbf_object(M = 11) # Parameters of the BPR model
opts$lambda        <- 1/20  # Regularisation term
opts$all_cells     <- length(region_dt)  # Total number of cells
opts$gpc_cov       <- 10      # Keep regions with enough GpC coverage
opts$BIC_thresh    <- 3       # Threshold in BIC difference between clusters
opts$tot_cl        <- 100     # maximum number of clusters
region_dt          <- lapply(region_dt, function(x) lapply(x, function(y){ 
  if(NROW(y) < opts$gpc_cov) return(NA) else return(y) }))

dt_N <- length(region_dt)
# Initialize data table
dt <- data.table(id=character(dt_N), clusters=numeric(dt_N), 
                 cell_clusters=list(), cell_profiles=list(), coverage = numeric(dt_N))

genes_len <- vector(mode = "numeric", length = opts$all_cells)
cell_names <- names(region_dt)
# Iterate over each cell
for (i in 1:dt_N){
  prev_BIC <- 1e100
  print(i)
  genes <- region_dt[[i]]             # Extract all genes for specific cell
  ind   <- which(is.na(genes))        # Filter to keep only GpC covered regions
  genes <- genes[-ind]
  
  uniq_genes <- length(genes)         # Number of unique genes having coverage
  cl_obj     <- list()
  for (k in 1:opts$tot_cl){  # Iterate over number of clusters
    # Cluster profiles
    cl_obj[[k]] <- bpr_cluster_wrap(x             = genes,
                                    K             = k,
                                    basis         = opts$basis_prof,
                                    lambda        = opts$lambda,
                                    em_max_iter   = 100,
                                    opt_itnmax    = 30,
                                    is_verbose    = TRUE,
                                    is_parallel   = FALSE,
                                    no_cores      = 3)
    # Check if we have a lower BIC
    if (prev_BIC - cl_obj[[k]]$BIC > opts$BIC_thresh){
      prev_BIC <- cl_obj[[k]]$BIC
    }else{
      # Extract cluster assignments of cells
      cell_clusters <- list(cell_clusters = cl_obj[[k-1]]$labels)
      # Extract cluster profiles
      cell_profiles <- list(cell_profiles = cl_obj[[k-1]]$w)
      set(dt, i, 1:5, as.list(c(cell_names[i], k-1, list(cell_clusters), list(cell_profiles), uniq_genes)))
      break
    }
    # In case every cell is a different cluster
    if (k == opts$tot_cl){
      cell_clusters <- list(cell_clusters = cl_obj[[k-1]]$labels)
      cell_profiles <- list(cell_profiles = cl_obj[[k-1]]$w)
      set(dt, i, 1:5, as.list(c(cell_names[i], k-1, list(cell_clusters), list(cell_profiles), uniq_genes)))
    }
  }
}
# Store them locally
saveRDS(dt, paste0(io$data_dir, "cell_wise_cluster_", annos, "_basis", opts$basis_prof$M, 
                   "_GpCcov", opts$gpc_cov, ".rds"))

# cols <- c("weight_dist", "sc_NLL", "bulk_NLL", "coverage", "cells")
# dt[, (cols) := lapply(.SD, as.numeric), .SDcols=cols]
# xs <- seq(-1, 1, length.out = 100)
# plot(x = xs, y = eval_probit_function(opts$basis_prof, xs, bulk_prof$w_opt[1:(opts$basis_prof$M + 1)]), 
#      type = "l", xlab = "region", ylab = "1 - Accessibility", ylim = c(0,1))
