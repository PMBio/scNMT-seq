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

model_sel <- function(region_dt, opts, cell_names){
  dt <- data.table(id=character(1), cl_prop=list(), cl_BIC=list(), coverage = numeric(1))
  genes <- region_dt                  # Extract all genes for specific cell
  ind   <- which(is.na(genes))        # Filter to keep only GpC covered regions
  genes <- genes[-ind]
  uniq_genes <- length(genes)         # Number of unique genes having coverage
  cl_obj = prop <- list()
  BIC <- vector(mode = "numeric", length = opts$tot_cl)
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
    
    # Extract cluster assignments of cells
    prop[[k]] <- list(cell_clusters = cl_obj[[k-1]]$pi_k)
    BIC[k] <- cl_obj[[k]]$BIC
  }
  set(dt, 1, 1:4, as.list(c(cell_names, list(prop), list(BIC), uniq_genes)))
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
opts$gpc_cov       <- 8         # Keep regions with enough GpC coverage
opts$BIC_thresh    <- 10        # Threshold in BIC difference between clusters
opts$tot_cl        <- 50        # Maximum number of clusters
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

# Parallel optimization
res <- parallel::mclapply(X = 1:dt_N,
                              FUN = function(y)
                                model_sel(region_dt = region_dt[[y]], opts = opts, cell_names = cell_names[y]),
                              mc.cores = no_cores)
# Stop parallel execution
parallel::stopCluster(cl)
doParallel::stopImplicitCluster()
  
dd <- rbindlist(res)
# Store them locally
saveRDS(dd, paste0(io$data_dir, "cell_wise_cl_allBIC_", annos, "_basis", opts$basis_prof$M,
                   "_GpCcov", opts$gpc_cov, ".rds"))
