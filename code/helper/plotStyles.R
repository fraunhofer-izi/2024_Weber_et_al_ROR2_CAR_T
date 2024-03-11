# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Colors
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
.inst = c("ggthemes", "scales") %in% installed.packages()
if (any(!.inst)) {
  install.packages(.cran_packages[!.inst], repos = "http://cran.rstudio.com/")
}

colorblind_vector <- colorRampPalette(c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF", "#0348A6"))
tableau.10 = ggthemes::tableau_color_pal("Tableau 10")(10)
tableau.20 = ggthemes::tableau_color_pal("Tableau 20")(20)
colors_stata =  ggthemes::stata_pal("s2color")(15)

# https://tradeblotter.wordpress.com/2013/02/28/the-paul-tol-21-color-salute/
tol2qualitative=c("#4477AA", "#CC6677")
tol3qualitative=c("#4477AA", "#DDCC77", "#CC6677")
tol4qualitative=c("#4477AA", "#117733", "#DDCC77", "#CC6677")
tol5qualitative=c("#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677")
tol6qualitative=c("#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677","#AA4499")
tol7qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#DDCC77", "#CC6677","#AA4499")
tol8qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677","#AA4499")
tol9qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499")
tol10qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#882255", "#AA4499")
tol11qualitative=c("#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#882255", "#AA4499")
tol12qualitative=c("#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#AA4466", "#882255", "#AA4499")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# THEME: GGPLOT
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
mytheme = function(base_size = 8, base_family = "") {
  half_line <- base_size/2
  theme_light(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "transparent",colour = NA),
      plot.background = element_rect(fill = "transparent",colour = NA),
      axis.ticks.length = unit(half_line / 2.2, "pt"),
      axis.ticks = element_line(colour = "black"),
      strip.background = element_rect(fill = NA, colour = NA),
      strip.text.x = element_text(size = rel(1), colour = "black"),
      strip.text.y = element_text(size = rel(1), colour = "black"),
      strip.text = element_text(size = rel(1), colour = "black"),
      axis.text = element_text(size = rel(1), colour = "black"),
      axis.title = element_text(size = rel(1), colour = "black"),
      legend.title = element_text(colour = "black", size = rel(1)),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = .3),
      legend.key.size = unit(1, "lines"),
      legend.text = element_text(size = rel(1), colour = "black"),
      legend.key = element_rect(colour = NA, fill = NA),
      legend.background = element_rect(colour = NA, fill = NA),
      plot.title = element_text(hjust = 0, face = "plain", colour = "black", size = rel(1)),
      plot.subtitle = element_text(colour = "black", size = rel(.85))
    )
}

mytheme_grid = function(base_size = 8, base_family = "") {
  half_line <- base_size/2
  theme_light(base_size = base_size, base_family = base_family) +
    theme(
      panel.background = element_rect(fill = "transparent",colour = NA),
      plot.background = element_rect(fill = "transparent",colour = NA),
      axis.ticks.length = unit(half_line / 2.2, "pt"),
      axis.ticks = element_line(colour = "black"),
      strip.background = element_rect(fill = NA, colour = NA),
      strip.text.x = element_text(size = rel(1), colour = "black"),
      strip.text.y = element_text(size = rel(1), colour = "black"),
      strip.text = element_text(size = rel(1), colour = "black"),
      axis.text = element_text(size = rel(1), colour = "black"),
      axis.title = element_text(size = rel(1), colour = "black"),
      legend.title = element_text(colour = "black", size = rel(1)),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = .3),
      legend.key.size = unit(1, "lines"),
      legend.text = element_text(size = rel(1), colour = "black"),
      legend.key = element_rect(colour = NA, fill = NA),
      legend.background = element_rect(colour = NA, fill = NA),
      plot.title = element_text(hjust = 0, face = "plain", colour = "black", size = rel(1)),
      plot.subtitle = element_text(colour = "black", size = rel(.85))
    )
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# MISC: GGPLOT
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
my.breaks <- function(my.vector, step) {
  my.min <- floor(min(my.vector))
  my.max <- ceiling(max(my.vector))
  my.seq <- seq(
    from = round(my.min / step) * step,
    to   = round(my.max / step) * step,
    by   = step
  )
  return(my.seq)
}

every_nth <- function(x,
                      nth,
                      empty = TRUE,
                      inverse = FALSE)
{
  if (!inverse) {
    if (empty) {
      x[1:nth == 1] <- ""
      x
    } else {
      x[1:nth != 1]
    }
  } else {
    if (empty) {
      x[1:nth != 1] <- ""
      x
    } else {
      x[1:nth == 1]
    }
  }
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# HEATMAP COLORS
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
hm.col = c(
  "#4477AA",
  "#66CCEE",
  "#228833",
  "#CCBB44",
  "#EE6677",
  "#AA3377",
  "#AE76A3" ,
  "#BBBBBB",
  "#000000"
)
degColors <- function(ann,
                      col_fun = FALSE,
                      con_values = c("grey80", "black"),
                      cat_values = c("orange", "steelblue"),
                      four_values = c("#44AA99", "#CC6677", "#DDCC77", "#4477AA"),
                      tol7qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#DDCC77", "#CC6677","#AA4499"),
                      palette = "Set2",
                      hm.col = hm.col) {
  col <- lapply(names(ann), function(a) {
    # if (class(ann[[a]][1]) == "numeric" | class(ann[[a]][1]) == "integer"){
    #     fn = colorRamp2(c(min(ann[[a]]),s
    #                       max(ann[[a]])),
    #                     con_values)
    #     if (col_fun){
    #         return(fn)
    #     }else{
    #         return(fn(ann[[a]]))
    #     }
    # }

    if (length(unique(ann[[a]])) < 3) {
      v <- cat_values[1:length(unique(ann[[a]]))]
      names(v) <- unique(ann[[a]])
      return(v)
    }
    if (length(unique(ann[[a]])) == 4) {
      v <- four_values[1:length(unique(ann[[a]]))]
      names(v) <- unique(ann[[a]])
      return(v)
    }
    if (length(unique(ann[[a]])) == 7) {
      v <- tol7qualitative[1:length(unique(ann[[a]]))]
      names(v) <- unique(ann[[a]])
      return(v)
    }
    if (length(unique(ann[[a]])) > 7 & length(unique(ann[[a]])) < 12) {
      v <- tol11qualitative[1:length(unique(ann[[a]]))]
      names(v) <- unique(ann[[a]])
      return(v)
    }
    if (length(unique(ann[[a]])) > 11) {
      v <- tol21rainbow[1:length(unique(ann[[a]]))]
      names(v) <- unique(ann[[a]])
      return(v)
    }
    v = hm.col[sort(unique(factor(ann[[a]])))]
    # v =  brewer.pal(length(unique(ann[[a]])), palette)
    names(v) <- sort(unique(ann[[a]]))
    v
  })
  names(col) <- names(ann)
  col
}
