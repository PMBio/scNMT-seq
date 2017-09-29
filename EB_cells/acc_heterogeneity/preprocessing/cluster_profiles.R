
# devtools::install_github("andreaskapou/BPRMeth-devel")
suppressPackageStartupMessages(library(BPRMeth))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(pracma))
suppressPackageStartupMessages(library(data.table))

# Initialize parameters
io            <- list(type = "acc")
io$base_dir   <- "scNMT_data_EB_dir"
io$data_dir   <- paste0(io$base_dir, "/", io$type, "/parsed/profiles/")

# Load saved objects for accessibility data
load(paste0(io$data_dir, io$type, "_400bp.RData"))
io$out_dir   <- paste0(io$base_dir, "/", io$type, "/parsed/profiles/")

annos              <- "prom_200_200"
opts$basis_prof    <- create_rbf_object(M = 13) # Parameters of the BPR model
opts$lambda        <- 1/20                      # Regularisation term
opts$all_cells     <- length(region_dt)         # Total number of cells
opts$cell_cov_prcg <- 0.6                       # Keep regions that are covered across many cells
opts$BIC_thresh    <- 2                         # Threshold in BIC difference between clusters
opts$gpc_cov       <- 10                        # Keep regions with enough GpC coverage
dt_N               <- NROW(anno_dt)             # Number of promoter regions
region_dt          <- lapply(region_dt, function(x) lapply(x, function(y){ 
  if(NROW(y) < opts$gpc_cov) return(NA) else return(y) }))
# Initialize data table
dt <- data.table(id=character(dt_N), cell_names=list(), anno=character(dt_N), 
                 clusters=numeric(dt_N), cell_clusters=list(), cell_profiles=list(), 
                 coverage = numeric(dt_N), cells = numeric(dt_N))

# Iterate over each genomic region
for (i in 1:dt_N){
  prev_BIC <- 1e100
  print(i)
  cells <- lapply(region_dt, "[[", i)  # Extract all cells for specific region
  ind   <- which(is.na(cells))         # Filter to keep only GpC covered regions
  cells <- cells[-ind]
  bulk  <- do.call(rbind, cells)       # Concatenate to obtain bulk data
  # If percentage of uncovered regions is less than threshold
  if (!is.null(bulk) & (length(ind) / opts$all_cells < opts$cell_cov_prcg)){
    uniq_cells <- length(cells)        # Number of unique cells having coverage
    cov        <- NROW(bulk)           # GpC coverage
    cell_names <- list(cell_names = names(cells))  # Cell names
    cl_obj     <- list()
    for (k in 1:length(cells)){        # Iterate over number of clusters
      # Cluster profiles
      cl_obj[[k]] <- bpr_cluster_wrap(x             = cells,
                                      K             = k,
                                      basis         = opts$basis_prof,
                                      lambda        = opts$lambda,
                                      em_max_iter   = 30,
                                      opt_itnmax    = 30,
                                      is_parallel   = FALSE)
      # Check if we have a lower BIC
      if (prev_BIC - cl_obj[[k]]$BIC > opts$BIC_thresh){
        prev_BIC <- cl_obj[[k]]$BIC
      }else{
        # Extract cluster assignments of cells
        cell_clusters <- list(cell_clusters = cl_obj[[k-1]]$labels)
        # Extract cluster profiles
        cell_profiles <- list(cell_profiles = cl_obj[[k-1]]$w)
        set(dt, i, 1:8, as.list(c(anno_dt$id[i], list(cell_names), annos, k-1, 
                                  list(cell_clusters), list(cell_profiles), cov, uniq_cells)))
        break
      }
      # In case every cell is a different cluster
      if (k == length(cells)){
        cell_clusters <- list(cell_clusters = cl_obj[[k-1]]$labels)
        cell_profiles <- list(cell_profiles = cl_obj[[k-1]]$w)
        set(dt, i, 1:8, as.list(c(anno_dt$id[i], list(cell_names), annos, k-1, 
                                  list(cell_clusters), list(cell_profiles), cov, uniq_cells)))
      }
    }
  }
}
# Store results
saveRDS(dt, paste0(io$out_dir, "cons_cluster_", annos, "_basis", opts$basis_prof$M, "_GpCcov", 
                   opts$gpc_cov, "_bic", opts$BIC_thresh,  "_cellcov", opts$cell_cov_prcg, ".rds"))
