LoadRNACounts = function(counts.file, name.slot = 4){
  # Load raw RNAseq counts from seqmonk report (tab delim text)
  d = fread(counts.file, colClasses = c("Chromosome" = "factor"), showProgress=F)
  n = colnames(d)[13:ncol(d)] %>%
    strsplit(split = "_") %>%
    map(~paste(.[name.slot], .[length(.)], sep = ".")) %>%
    unlist()
  setnames(d, 13:ncol(d), n)
  d = d[, Gene := make.names(Probe, unique = TRUE)] # make sure each gene name is unique
  return(d)
}
