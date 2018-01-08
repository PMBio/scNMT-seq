

FindFiles = function(dir.list, cg.pattern = "CpG_report.gz$", gc.pattern = "GpC_report",
                     name.slot = c(13, 22, 23), read.slot = 15){
  cg = unlist(llply(dir.list, function(x) paste0(x, list.files(x, pattern = cg.pattern))))
  gc = unlist(llply(dir.list, function(x) paste0(x, list.files(x, pattern = gc.pattern))))
  
  cg = unlist(llply(cg, function(x) ifelse(isGzipped(x), paste("zcat", x), x)))
  gc = unlist(llply(gc, function(x) ifelse(isGzipped(x), paste("zcat", x), x)))
  
  files = c(cg, gc)
  names = unlist(llply(files, function(x) paste(unlist(strsplit(x, split = "_"))[name.slot], collapse = ".")))
  read = unlist(llply(files, function(x) unlist(strsplit(x, split = "_"))[read.slot]))
  context = unlist(llply(names, function(x) ifelse(any(grep("GpC", x)), "GC", "CG")))
  type = paste(context, read, sep = ".")
  names = unlist(llply(names, function(x) paste(unlist(strsplit(x, split = "\\."))[c(1, 4)], collapse = ".")))
  files = data.table(name = names, file = files, context = context, read = read )
  files = dcast(files, name ~ type, value.var = "file")
  files$n = row.names(files)
  return(files)
}

ReadGCandCG = function(gc.R1.file, gc.R2.file, cg.R1.file, cg.R2.file, name){
  gc = list(
    r1 = fread(gc.R1.file, select = c(1:2, 4:6), colClasses = list(factor = 1L), showProgress = FALSE),
    r2 = fread(gc.R2.file, select = c(1:2, 4:6), colClasses = list(factor = 1L), showProgress = FALSE))
  gc = rbindlist(gc, use.names = TRUE)[
    V4 > 0 | V5 > 0]
  setkey(gc, V6)
  gcg = gc["CG"]
  gcg[, c("V4", "V5", "V6") := NULL][, GC := "GCG"]
  gc = gc[!"CG"]
  setkey(gc, V1)
  gc = gc[!"MT"][ # remove mitochondria
          , V6 := NULL][
          , lapply(.SD, sum), .(V1, V2)][ # combine read 1 and read 2
          , V6 := 100 * V4 / (V4 + V5)][ # calculate mean methylation
          , c("V4", "V5") := NULL]
  setnames(gc, 1:3, c("chr", "bp", name))
  
  cg = list(
    r1 = fread(cg.R1.file, select = c(1:2, 4:6), colClasses = list(factor = 1L), showProgress = FALSE),
    r2 = fread(cg.R2.file, select = c(1:2, 4:6), colClasses = list(factor = 1L), showProgress = FALSE))
  cg = rbindlist(cg, use.names = TRUE)[
    V4 > 0 | V5 > 0]
  setkey(cg, V6)
  cg = cg["CG"]
  setkey(cg, V1, V2)
  setkey(gcg, V1, V2)
  cg = merge(cg, gcg, all = TRUE)
  setkey(cg, GC)
  cg = cg[!"GCG"][, GC := NULL]
  
  setkey(cg, V1)
  cg = cg[!"MT"]
  cg[, V6 := NULL]
  
  cg = cg[, lapply(.SD, sum), .(V1, V2)][
          , V6 := 100 * V4 / (V4 + V5)][
          , c("V4", "V5") := NULL]
  
  setnames(cg, 1:3, c("chr", "bp", name))
  
  return(list(gc = gc, cg = cg))
  
}

ReadCG = function(cg.R1.file, cg.R2.file, name){
  cg = list(
    r1 = fread(cg.R1.file, select = c(1:2, 5:6), colClasses = list(factor = 1L), showProgress = FALSE),
    r2 = fread(cg.R2.file, select = c(1:2, 5:6), colClasses = list(factor = 1L), showProgress = FALSE))
  cg = rbindlist(cg, use.names = TRUE)
  
  
  setkey(cg, V1)
  cg = cg[!"MT"]
  
  # combine Read 1 and read 2 then compute mean methylation
  setkey(cg, V1, V2)
  cg = cg[, lapply(.SD, sum), .(V1, V2)]
  cg[, V7 := 100*V6/(V5+V6)]
  cg[, c("V5", "V6") := NULL]
  
  
  setnames(cg, 1:3, c("chr", "bp", name))

  return(cg)
  
}


QCRawMeth = function(cg, gc, cell){
  qc = list(cg[, .(CpGs = .N, 
                   cg.binary = cg[get(cell) ==100 | get(cell) ==0, .N],
                   cg.mean = mean(get(cell), na.rm = TRUE))],
            gc[, .(GpCs = .N, 
                   gc.binary = gc[get(cell) ==100 | get(cell) ==0, .N],
                   gc.mean = mean(get(cell), na.rm = TRUE))])
  qc = cbind(qc[[1]], qc[[2]])
  qc[, cell := cell]
  
  return(qc)
}

Quantitate = function(overlap.data){
  # Quantitate methylation over regions
  # Args: overlap.data = data.table of methylation calls subseted within regions
  # Returns: list of 2 data.table - meth = mean methylation, n = number of CpGs covered
  names = colnames(overlap.data)
  cell = names[!names %in% c("chr", "Start", "End", "Probe", "TSS", "bp", "Strand")]
  setkey(overlap.data, bp, Start, End)
  overlap.data = overlap.data[bp >= Start & bp <= End]
  overlap.data[, bp  := NULL]
  if(any(names %in% "Probe")) {
    meth = overlap.data[, lapply(.SD, mean, na.rm = TRUE), .(chr, Start, End, Probe), .SDcol = cell]
    n = overlap.data[, lapply(.SD, function(x) sum(!is.na(x))), .(chr, Start, End, Probe), .SDcol = cell]
  } else {
    meth = overlap.data[, lapply(.SD, mean, na.rm = TRUE), .(chr, Start, End), .SDcol = cell]
    n = overlap.data[, lapply(.SD, function(x) sum(!is.na(x))), .(chr, Start, End), .SDcol = cell]
  }
  return(list(meth = meth, n = n))
  
}

MergeAll = function(x, y) merge(x, y, all = TRUE)

MergeCells = function(data){
  # Merge cells / samples from list format into two data.tables of meth and n
  meth = llply(data, function(x) x$meth)
  cols = colnames(meth[[1]])
  cols = cols[cols %in% c("chr", "Start", "End", "Probe")]
  l_ply(meth, setkeyv, cols)
  n = llply(data, function(x) x$n)
  l_ply(n, setkeyv, cols)
  data = list(meth = Reduce("MergeAll", meth),
              n = Reduce("MergeAll", n))
  return(data)
}


convert_chr_format <- function(chr, to) {
  # Function to convert the chr from short to long format and viceversa
  # to: "short" or "long"
  chr <- as.character(chr)
  stopifnot(to %in% c("short","long"))
  short_alphabet <- c(1:19,"X","Y","MT")
  long_alphabet <- paste("chr",short_alphabet,sep="")
  if (to == "short") {
    if (all(chr %in% short_alphabet)) { 
      return(chr) 
    } else {
      stopifnot(all(chr %in% long_alphabet))
      names(short_alphabet) <- long_alphabet
      return(unname(short_alphabet[chr]))
    }
  }
  if (to == "long") {
    if (all(chr %in% long_alphabet)) { 
      return(chr) 
    } else {
      stopifnot(all(chr %in% short_alphabet))
      names(long_alphabet) <- short_alphabet
      return(unname(long_alphabet[chr]))
    }
  }
}


extract_seq <- function(chr, pos1, pos2, strand="+") {
  chr <- convert_chr_format(chr, to="long")
  d <- getSeq(Mmusculus, chr, pos1, pos2) # Required BSgenome package
  d[strand=="-"] <- complement(d[strand=="-"])
  return(d)
}

