theme_pub <- function() {
    theme(
    plot.title = element_text(size=20, hjust=0.5, margin=margin(0,0,20,0)),
    axis.title.y = element_text(colour="black", size=10, vjust=1.5),
    axis.title.x = element_text(colour="black", size=10, vjust=1.5, margin=margin(15,0,0,0)),
    axis.text.x = element_text(colour="black",size=rel(1.2)),
    axis.text.y = element_text(colour="black",size=rel(1.2)),
    # axis.line = element_line(colour="black", size=rel(0.7)),
    # axis.ticks.x = element_line(colour="black", size=rel(0.8)),
    # axis.ticks.y = element_blank(),
    legend.position="right",
    legend.title=element_blank(),
    legend.key       = element_rect(fill = "white", colour = "white"),
    legend.background = element_rect(fill="white", color="white"),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_line(color="grey")
  )
}

theme_bw <- function() {
    theme(
      # white background and dark border
      panel.background = element_rect(fill = "white", colour = NA),
      panel.border     = element_rect(fill = NA, colour = "grey20"),
      # make gridlines dark, same contrast with white as in theme_grey
      panel.grid.major = element_line(colour = "grey92"),
      panel.grid.minor = element_line(colour = "grey92", size = rel(0.5)),
      # contour strips to match panel contour
      strip.background = element_rect(fill = "grey85", colour = "grey20"),
      # match legend key to background
      legend.key       = element_rect(fill = "white", colour = NA),

      complete = TRUE
    )
}


sample_colors <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}