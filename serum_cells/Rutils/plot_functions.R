
PlotMeth = function(data, sample, title, contexts = c("GpC", "CpG"), stratify = "context", tile = 30, xlim = 1000){
  # plot methylation data over given region
  # Args: data = data.table with methylation values
  #       sample = name of sample to plot - must be a column name in data
  #       title = title of plot
  #       contexts = contexts to plot - either "CpG", "GpC" or c("CpG", "GpC")
  #       tile = window to average methylation values over
  setkey(data, context)
  to.plot = copy(data)[context %in% contexts][, dist := tile * round(dist / tile)][
    , lapply(.SD, mean, na.rm = TRUE), .(dist, get(stratify)), .SDcol = sample] %>%
    setnames(sample, "Methylation") %>%
    setnames("get", stratify)
  p = ggplot(to.plot, aes(dist, Methylation, colour = get(stratify))) + 
    geom_line() + 
    ggtitle(title) + 
    coord_cartesian(xlim = c(-xlim, xlim)) +
    xlim(-xlim, xlim) +
    theme(legend.title=element_blank()) +
    xlab("Genomic distance") +
    scale_colour_discrete(labels = c("CpG Methylation",
                                     "GpC Accessibility"))
  return(p)
}

