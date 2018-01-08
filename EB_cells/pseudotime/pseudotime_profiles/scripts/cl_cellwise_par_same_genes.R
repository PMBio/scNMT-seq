# ------------------------------------------
# Set working directory and load libraries
# ------------------------------------------
if (interactive()){
  cur.dir <- dirname(parent.frame(2)$ofile)
  setwd(cur.dir)
}

set.seed(1)
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(BPRMeth))
suppressPackageStartupMessages(library(pracma))
suppressPackageStartupMessages(library(data.table))

model_sel <- function(region_dt, opts, cell_names){
  dt <- data.table(id=character(1), clusters=numeric(1), 
                   cell_clusters=list(), cell_profiles=list(), coverage = numeric(1))
  prev_BIC <- 1e100
  genes <- region_dt                  # Extract all genes for specific cell
  ind   <- which(is.na(genes))        # Filter to keep only GpC covered regions
  genes <- genes[-ind]
  gene_samples <- sample(length(genes), opts$gene_num)
  # genes <- genes[1:opts$gene_num]
  genes <- genes[gene_samples]
  uniq_genes <- length(genes)         # Number of unique genes having coverage
  cl_obj     <- list()
  for (k in 5:opts$tot_cl){  # Iterate over number of clusters
    # Cluster profiles
    cl_obj[[k]] <- bpr_cluster_wrap(x             = genes,
                                    K             = k,
                                    basis         = opts$basis_prof,
                                    lambda        = opts$lambda,
                                    em_max_iter   = 150,
                                    opt_itnmax    = 50,
                                    epsilon_conv  = 1e-03,
                                    is_verbose    = TRUE,
                                    is_parallel   = FALSE)
    # Check if we have a lower BIC
    if (prev_BIC - cl_obj[[k]]$BIC > opts$BIC_thresh){
      prev_BIC <- cl_obj[[k]]$BIC
    }else{
      # Extract cluster assignments of cells
      cell_clusters <- list(cell_clusters = cl_obj[[k-1]]$labels)
      # Extract cluster profiles
      cell_profiles <- list(cell_profiles = cl_obj[[k-1]]$w)
      set(dt, 1, 1:5, as.list(c(cell_names, k-1, list(cell_clusters), list(cell_profiles), uniq_genes)))
      break
    }
    # In case every cell is a different cluster
    if (k == opts$tot_cl){
      cell_clusters <- list(cell_clusters = cl_obj[[k]]$labels)
      cell_profiles <- list(cell_profiles = cl_obj[[k]]$w)
      set(dt, 1, 1:5, as.list(c(cell_names, k, list(cell_clusters), list(cell_profiles), uniq_genes)))
    }
  }
  return(dt)
}


# Initialize parameters
io            <- list(type = "acc")
io$base_dir   <- "../datasets/nmt-seq-data-eb"
io$data_dir   <- paste0(io$base_dir, "/", io$type, "/parsed/profiles/")
# Load saved objects for accessibility data
load(paste0(io$data_dir, io$type, "_400bp.RData"))
io$data_dir   <- paste0(io$base_dir, "/", io$type, "/parsed/profiles/")

dt_N               <- length(region_dt)   # Number of cells
annos              <- "prom_200_200"
no_cores           <- dt_N
opts$basis_prof    <- create_rbf_object(M = 11) # Parameters of the BPR model
opts$lambda        <- 1/20  # Regularisation term
opts$all_cells     <- length(region_dt)  # Total number of cells
opts$gpc_cov       <- 8        # Keep regions with enough GpC coverage
opts$BIC_thresh    <- 10       # Threshold in BIC difference between clusters
opts$tot_cl        <- 100      # Maximum number of clusters
region_dt          <- lapply(region_dt, function(x) lapply(x, function(y){ 
  if(NROW(y) < opts$gpc_cov) return(NA) else return(y) }))

cell_names <- names(region_dt)

# If number of cores is not given
if (is.null(no_cores)){ no_cores <- parallel::detectCores() - 2
}else{ if (no_cores >= parallel::detectCores()){ no_cores <- parallel::detectCores() - 1 } }
if (is.na(no_cores)){ no_cores <- 2 }
# Create cluster object
cl <- parallel::makeCluster(no_cores)
doParallel::registerDoParallel(cl)

# Get how many genes are at most in each cell
tot_genes <- c()
for (i in 1:dt_N){
  genes <- region_dt[[i]]             # Extract all genes for specific cell
  ind   <- which(is.na(genes))        # Filter to keep only GpC covered regions
  genes <- genes[-ind]
  tot_genes <- c(tot_genes, length(genes))
}
opts$gene_num <- min(tot_genes)
# Parallel optimization
res <- parallel::mclapply(X = 1:dt_N,
                              FUN = function(y)
                                model_sel(region_dt = region_dt[[y]], opts = opts, cell_names = cell_names[y]),
                              mc.cores = no_cores)
# Stop parallel execution
parallel::stopCluster(cl)
doParallel::stopImplicitCluster()
  
dd <- rbindlist(res) # Get a data.table
# Store them locally
saveRDS(dd, paste0(io$data_dir, "cell_wise_cl_sameGenes_", annos, "_basis", opts$basis_prof$M,
                   "_GpCcov", opts$gpc_cov, ".rds"))
