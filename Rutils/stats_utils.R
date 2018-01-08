theme_barplot_pub <- function() {
  p <- theme(
    legend.position="none",
    axis.title.y = element_text(colour="black", size=20),
    axis.text.x = element_text(colour="black", angle=90, size=10, vjust=0.5, hjust=1.0),
    axis.ticks.x = element_line(colour="black")
  )
  return(p)
}

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

barplot_theme <- function() {
  p <- theme(
    plot.title = element_text(size=30, hjust=0.5),
    # axis.title.x = element_text(colour="black", size=25, vjust=1.5),
    axis.title.x = element_blank(),
    axis.title.y = element_text(colour="black", size=20),
    # axis.text.x = element_text(colour="black",size=rel(1.6)),
    axis.text.y = element_text(colour="black",size=rel(1.5)),
    axis.line = element_line(colour="black", size=rel(0.7)),
    axis.ticks.x = element_line(colour="black", size=rel(1.0)),
    axis.ticks.y = element_line(colour="black", size=rel(1.0)),
    legend.position="none",
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank()
  )
}

gg_qqplot = function(xs, perm_xs, ci=0.95, title = "Quantile-quantile plot of p-values") {
  N = length(xs)
  df = data.frame(observed=-log10(sort(xs)),
                  permuted=-log10(sort(perm_xs)),
                  expected=-log10(1:N / N),
                  cupper=-log10(qbeta(ci,     1:N, N - 1:N + 1)),
                  clower=-log10(qbeta(1 - ci, 1:N, N - 1:N + 1)))
  
  log10Pe = expression(paste("Expected -log"[10], plain(P)))
  log10Po = expression(paste("Observed -log"[10], plain(P)))
  ggplot(df) +
    geom_point(aes(expected, observed), shape=1, size=2.5) +
    geom_point(aes(expected, permuted), shape=2, size=1.5, color = "coral") +
    geom_abline(intercept=0, slope=1, alpha=0.5, color = "red") +
    geom_line(aes(expected, cupper), linetype=2) +
    geom_line(aes(expected, clower), linetype=2) +
    xlab(log10Pe) +
    ylab(log10Po) + 
    labs(title=title) +
    theme(
      plot.title=element_text(size=rel(1.7), hjust=0.5),
      plot.subtitle=element_text(size=rel(1.6), face="italic", hjust=0.5),
      axis.text=element_text(size=rel(1.2), color='black'),
      axis.title=element_text(size=rel(1.3), color='black'),
      panel.border=element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank()
    )
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

impute <- function(d, margin) {
  if (margin == 1)
    means <- rowMeans(d, na.rm=T)
  else if (margin == 2)
    means <- colMeans(d, na.rm=T)
  else
    stop("Margin has to be either 1 (rows) or 2 (cols)")
  
  if (any(is.na(means))) {
    stop('Insufficient data for mean imputation!')
  }
  
  for (i in 1:length(means)) {
    if (margin == 1)
      d[i,is.na(d[i,])] <- means[i]
    else if (margin == 2)
      d[is.na(d[,i]), i] <- means[i]
  }
  return (d)
}

pca <- function(d, center=T, scale=F) {
  # center and scale the data
  d <- scale(d, center=center, scale=scale)
  d <- t(d)
  # perform eigenalue decomposition on the cov matrix via svd on the data
  # (it is the same)
  s <- svd(d)
  # eigenvalues (principal components)
  vec <- s$u
  rownames(vec) <- rownames(d)
  # square the eigenvalues to obtain variance explained
  val <- s$d**2
  # calculate percentage of variation for each PC
  val <- val / sum(val)
  return (list(vec=vec, val=val))
}

# Function to plot the projected data onto PCX and PCY
plot_pca_vec <- function(pc_vec, x=1, y=2, title="") {
  t <- data.frame(sample=factor(rownames(pc_vec)), pcx=pc_vec[,x], pcy=pc_vec[,y])
  cols <- samples_colors(as.vector(t$sample))
  p <- ggplot(t, aes(x=pcx, y=pcy)) +
    geom_point(aes(color=sample), size=2) +
    scale_color_manual(values=cols) +
    ggtitle(title) +
    #geom_text(aes(label=sample), vjust=-.4, hjust= .3, size=3) +
    xlab(sprintf('pc%d', x)) + ylab(sprintf('pc%d', y)) +
    guides(color=F)# + theme_pub()
  return (p)
}

# Scree plot
plot_pca_val <- function(pc_val, title="") {
  t <- data.frame(pc=1:length(pc_val), val=pc_val)
  p <- ggplot(t, aes(x=pc, y=val)) +
    geom_bar(stat='identity', fill='salmon', color='black') +
    ggtitle(title) +
    xlab('principle component') +
    ylab('% variance explained')# + theme_pub()
  return (p)
}


barplot_theme <- function() {
  p <- theme(
    plot.title = element_text(size=30, hjust=0.5),
    axis.title.y = element_text(colour="black", size=25, vjust=1.5),
    axis.title.x = element_text(colour="black", size=25, vjust=1.5),
    axis.text.x = element_text(colour="black",size=rel(1.6)),
    axis.text.y = element_text(colour="black",size=rel(1.6)),
    axis.line = element_line(colour="black", size=rel(0.7)),
    axis.ticks.x = element_line(colour="black", size=rel(1.0)),
    axis.ticks.y = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank()
  )
}

barPlot <- function(data, xlabel="", ylabel="") {
  if (ncol(data) == 3) {
    p <- ggplot(data, mapping=aes_string(x=names(data)[1], y=names(data)[2], fill=names(data)[3]))
  } else if (ncol(data)==2) {
    p <- ggplot(data, mapping=aes_string(x=names(data)[1], y=names(data)[2]))
  }
  p + ggtitle("") + geom_bar(stat='identity', position="dodge") +
    xlab(xlabel) + ylab(ylabel) + barplot_theme()
}



boxplot_theme <- function() {
  p <- theme(
    plot.title = element_text(size=30, hjust=0.5, margin=margin(0,0,20,0)),
    axis.title.y = element_text(colour="black", size=25, vjust=1.5),
    axis.title.x = element_text(colour="black", size=25, vjust=1.5),
    axis.text.x = element_text(colour="black",size=rel(1.7)),
    axis.text.y = element_text(colour="black",size=rel(1.7)),
    axis.line = element_line(colour="black", size=rel(0.7)),
    axis.ticks.x = element_line(colour="black", size=rel(0.8)),
    axis.ticks.y = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    legend.position="top",
    legend.text=element_text(size=15),
    legend.title=element_blank(),
    legend.background=element_blank(),
    panel.border = element_blank()
  )
}

boxPlot2 <- function(data, xlabel="", ylabel="") {
  p <- ggplot(data, mapping=aes_string(x=names(data)[1], y=names(data)[2], fill=names(data)[3], group=names(data)[3])) + 
    ggtitle("") +
    geom_boxplot(alpha=0.8, coef=0, outlier.shape=NA) +
    xlab(xlabel) + ylab(ylabel) +
    boxplot_theme()
  
  return(p)
}

boxPlot <- function(data, xlabel="", ylabel="") {
  # p <- ggplot(data, mapping=aes_string(x=names(data)[1], y=names(data)[2], fill=names(data)[1])) + 
  p <- ggplot(data, mapping=aes_string(x=names(data)[1], y=names(data)[2])) + 
    ggtitle("") +
    geom_boxplot(alpha=0.8, coef=0, outlier.shape=NA) +
    xlab(xlabel) + ylab(ylabel) +
    coord_flip() +
    boxplot_theme()
  return(p)
}


stats_tbldf <- function(d) {
  s <- summarise(d,
                 mean=mean(rate, na.rm=T),
                 var=var(rate, na.rm=T),
                 wtd_mean=weighted.mean(rate, weight, na.rm=T),
                 wtd_var=wtd.var(rate, weight, na.rm=T),
                 weight=sum(weight)
  )
  return (s)
}

stats_datatable <- function(d, colsby) {
  s <- d[,.(mean=mean(rate, na.rm=T), wtd_mean=weighted.mean(rate, weight, na.rm=T), 
       var=var(rate, na.rm=T), wtd_var=wtd.var(rate, weight, na.rm=T), 
       weight=sum(weight) ), by=colsby]
  return (s)
}

scatterPlot <- function(data, xlabel="", ylabel="", title="", subtitle="", lm=FALSE) {
  p <- ggplot(tmp, aes_string(x=names(data)[1], y=names(data)[2])) +
    labs(x=xlabel, y=ylabel, title=title, subtitle=subtitle) +
    # ggtitle(title) + xlab(xlab) + ylab(ylab) +
    geom_point() +
    theme(
      plot.title=element_text(size=rel(2.0), hjust=0.5),
      plot.subtitle=element_text(size=rel(1.6), face="italic", hjust=0.5),
      axis.text=element_text(size=rel(1.3), color='black'),
      axis.title=element_text(size=rel(1.5), color='black'),
      panel.border=element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank()
    )
  
  if (lm)
    p <- p + geom_smooth(method="lm")
  
  print(p)
}

corPlot <- function(data, title="") {
  # columns: r, p, sig
  p <- ggplot(data, aes(x=r, y=p)) + 
    ggtitle(title) +
    geom_vline(xintercept=0, colour="orange") +
    geom_vline(aes(xintercept = mean(r)), colour="blue", linetype="dashed") +
    geom_point(aes(color=sig)) +
    scale_color_manual(values=c("black","red")) +
    xlab("Pearson Correlation coefficient") + ylab("-log10 p-value") +
    theme(
      plot.title=element_text(size=15, face='bold', margin=margin(0,20,0,0), hjust=0.5),
      legend.position="none",
      axis.text=element_text(size=rel(1.2), color='black'),
      axis.title=element_text(size=rel(1.2), color='black'),
      panel.border=element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank()
    )
  return(p)
}

personalised_theme <- function() {
  theme(
    axis.text.x = element_text(colour="black",size=rel(1.7)),
    axis.text.y = element_text(colour="black",size=rel(1.5), vjust=0.5),
    axis.title.x = element_text(colour="black",size=15),
    axis.ticks.y = element_line(colour="black", size=rel(1.2)),
    legend.position="none"
  )
}
