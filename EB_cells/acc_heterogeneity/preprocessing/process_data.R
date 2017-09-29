
# devtools::install_github("andreaskapou/BPRMeth-devel")
suppressPackageStartupMessages(library(BPRMeth))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(pracma))
suppressPackageStartupMessages(library(data.table))

# Input/ Output
io            <- list(type = "acc")
io$base_dir   <- "scNMT_data_EB_dir"
io$annos_dir  <- paste0(io$base_dir, "/features/filt")
io$meta_file  <- paste0(io$base_dir, "/sample_sheet.csv")
io$rna_file   <- paste0(io$base_dir, "/rna/parsed/sceset.rds")
io$data_dir   <- paste0(io$base_dir, "/", io$type, "/raw/merged/binarised")

opts <- list()
opts$annos  <- c("prom_200_200")
opts$upstream <- -200
opts$downstream <- 200

# Load cell metadata
metadata <- fread(io$meta_file) %>% .[,c("sample", "culture", "pass_metQC", "pass_accQC", "pass_rnaQC")] %>%
  na.omit() %>% .[pass_metQC == TRUE & pass_accQC ==  TRUE & pass_rnaQC == TRUE]

# Read scRNA-Seq data as SCESet object
sceset <- readRDS(file = io$rna_file)
# Create data.table
rna <- exprs(sceset) %>% t %>% as.data.table(keep.rownames = "sample") %>% melt(id.vars = "sample", value.name = "expr", variable.name = "gene")
# Extract gene coordinates metadata to do the overlap
rna_metadata <- fData(sceset) %>% tibble::rownames_to_column("gene") %>% as.data.table %>% .[,c("chr", "start", "end", "gene", "ens_id")]  %>% 
  .[, chr := as.factor(sub("chr", "", chr))] %>% setnames("ens_id", "id")

# Read annotation data
anno_dt <- lapply(opts$annos, function(anno) fread(sprintf("%s/%s.bed", io$annos_dir, anno), 
  colClasses = c("character", "integer", "integer", "character", "character", "character"))[, c(1, 2, 3, 4, 5, 6)]) %>%  
  rbindlist %>% setnames(c("chr", "start", "end", "strand", "id", "anno")) %>% 
  .[,c("anno", "chr") := list(as.factor(anno), as.factor(sub("chr","",chr)))] %>% subset(id %in% rna_metadata$id) 

# Extract the annotation names
opts$annos_names <- unique(opts$annos)
# Extract unique sample names
opts$sample_names <- unique(metadata$sample)
# Create regions
anno_dt <- anno_dt %>% .[, c("tss") := list(floor((start + end)/2))] %>% 
  .[, c("start", "end") := list(tss + opts$upstream, tss + opts$downstream)] %>% setkey(anno, chr, start, end)

# Create genomic regions for each annotation
for (ann in opts$annos_names){
  region_dt <- list()
  # Iterate over each cell
  for (cell in opts$sample_names){
    bs_data <- fread(input = sprintf("zcat < %s/%s.tsv.gz", io$data_dir, cell), sep = "\t", header = TRUE, 
                     stringsAsFactors = TRUE, showProgress = FALSE) %>%
      setnames(c("chr", "pos", "rate")) %>%
      .[, c("start", "end", "meth_reads", "pos") := list(pos, pos, ifelse(rate > 50, 1, 0), NULL)] %>%
      .[,c("chr", "start", "end", "meth_reads")] %>% setkey(chr, start, end)
    bs_data <- GRanges(bs_data)
    
    # Create promoter methylation regions
    region_dt[[cell]] <- create_methyl_region(bs_data = GRanges(bs_data),
                                              prom_region = GenomicRanges::GRanges(anno_dt[ann]),
                                              cpg_density = 3,
                                              sd_thresh = -1,
                                              ignore_strand = TRUE,
                                              is_single_cell = TRUE)
  }
  save(region_dt, anno_dt, opts, io, file = paste0(io$base_dir, "/", io$type,"/parsed/profiles/", io$type, "_400bp.RData"))
}
