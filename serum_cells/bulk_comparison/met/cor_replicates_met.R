io <- list()
io$met.data1 <- "/Users/ricard/data/Ficz_2013/BS/rep1_parsed_mm10.txt.gz"
io$met.data2 <- "/Users/ricard/data/Ficz_2013/BS/rep2_parsed_mm10.txt.gz"
io$met.data3 <- "/Users/ricard/data/Ficz_2013/BS/rep3_parsed_mm10.txt.gz"

met <- list(
  rep1 = fread(paste0("zcat < ",io$met.data1)),
  rep2 = fread(paste0("zcat < ",io$met.data2)),
  rep3 = fread(paste0("zcat < ",io$met.data3))
)

r <- matrix(nr=length(met), nc=length(met)); diag(r) <- 1
# n <- matrix(nr=length(opts$cells), nc=length(opts$cells)); diag(r) <- 1
rownames(r) <- names(met); colnames(r) <- names(met)
for (i in 1:length(met)) {
  data_i <- met[[i]] %>% .[rate%in%c(0,100)]
  for (j in i:length(met)) {
    if (i!=j) {
      data_j <- met[[j]] %>% .[rate%in%c(0,100)]
      data <- merge(data_i,data_j, by=c("chr","pos"))
      r[i,j] <- r[j,i] <- mean(data$rate.x==data$rate.y)
      # n[i,j] <- n[j,i] <- nrow(data)
    }
  }
}
