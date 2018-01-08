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
load(paste0(io$data_dir, io$type, "_200bp.RData"))
io$data_dir   <- paste0(io$base_dir, "/", io$type, "/parsed/profiles/")

annos              <- "prom_100_100"
opts$basis_prof    <- create_rbf_object(M = 9)  # Parameters of the BPR model
opts$lambda        <- 1/20  # Regularisation term
opts$all_cells     <- length(region_dt)  # Total number of cells
opts$cell_cov_prcg <- 0.5   # Keep regions that are covered across many cells
opts$gpc_cov       <- 10    # Keep regions with enough GpC coverage
region_dt          <- lapply(region_dt, function(x) lapply(x, function(y){ 
  if(NROW(y) < opts$gpc_cov) return(NA) else return(y) }))
dt_N               <- NROW(anno_dt) # Number of promoter regions
# Initialize data table
dt <- data.table(id=character(dt_N), cell_names=list(), anno=character(dt_N), weight_dist=numeric(dt_N), 
                 sc_NLL=numeric(dt_N), bulk_NLL = numeric(dt_N), coverage = numeric(dt_N), cells = numeric(dt_N))

# Iterate over each region
for (i in 1:dt_N){
  print(i)
  cells <- lapply(region_dt, "[[", i)  # Extract all cells for specific region
  ind  <- which(is.na(cells))    # Filter to keep only GpC covered regions
  cells <- cells[-ind]
  bulk <- do.call(rbind, cells)  # Concatenate to obtain bulk data
  # If percentage of uncovered regions is less than threshold
  if (!is.null(bulk) & (length(ind) / opts$all_cells < opts$cell_cov_prcg)){
    # Learn bulk profile
    bulk_prof <- bpr_optim(x             = bulk,
                           basis         = opts$basis_prof,
                           lambda        = opts$lambda,
                           fit_feature   = "NLL",
                           cpg_dens_feat = TRUE,
                           opt_itnmax    = 50)
    # Learn profile for each single cell
    sc_prof <- bpr_optim(x             = cells,
                         basis         = opts$basis_prof,
                         lambda        = opts$lambda,
                         fit_feature   = "NLL",
                         cpg_dens_feat = FALSE,
                         opt_itnmax    = 50, 
                         is_parallel   = FALSE)
    
    sc_NLL <- sum(sc_prof$W_opt[, opts$basis_prof$M + 2])  # NLL of single cell data
    bulk_NLL <- bulk_prof$w_opt[opts$basis_prof$M + 2]     # NLL of bulk data 
    cov <- bulk_prof$w_opt[opts$basis_prof$M + 3]          # GpC coverage of bulk data
    uniq_cells <- length(cells)    # Number of unique cells having coverage
    # Euclidean distance of BPR weights
    weight_dist <- sum(dist(pnorm(sc_prof$W_opt[, 1:opts$basis_prof$M + 1]), method = "euclidean"))
    cell_names <- list(cell_names = names(cells))          # Cell names
    # Update row of data.table
    set(dt, i, 1:8, as.list(c(anno_dt$id[i], list(cell_names), annos, round(weight_dist, 2), 
                              round(sc_NLL, 2), round(bulk_NLL, 2), cov, uniq_cells)))
  }
}
# Store them locally
saveRDS(dt, paste0(io$data_dir, "cons_prof_", annos, "_basis", opts$basis_prof$M, "_GpCcov", 
                   opts$gpc_cov, "cellcov", opts$cell_cov_prcg, ".rds"))
