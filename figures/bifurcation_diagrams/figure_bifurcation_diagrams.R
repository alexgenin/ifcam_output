# 
# 
# This file contains code that generate bifurcation diagrams for all models. 
# 
# It generates 2D planes with contours for the upper branch (see Flo's 
#   manuscript, Fig 2.) and cross-sections chose for certain thresholds.
# 

# Some prerequisite before anything else
# Set proxy if needed

# Load libraries
library(ggplot2)
library(grid) 
library(gridExtra) # for grid.arrange

# Define output path 
# /!\ NB: folders must have a trailing /
working_directory  <- '/home/alex/work/2014-2015/SpatialStress/ifcam_output/figures/'
output_figure_path <- "/home/alex/work/2014-2015/SpatialStress/ifcam_output/figures/bifurcation_diagrams/" # same folder

# Define data paths 
# /!\ NB: folders must have a trailing /
data_folder <- '/home/alex/work/2014-2015/SpatialStress/ifcam_output/data/'
# Use this folder if you do not have the results on your computer: it will fetch
#   the latest ones from github.
# data_folder <- "https://github.com/alexgenin/ifcam_output/blob/master/data/" 
files <- list(musselbed = paste0(data_folder, "result_musselbed_processed.rda"),
              grazing   = paste0(data_folder, "result_grazing_processed.rda"),
              forestgap = paste0(data_folder, "result_forestgap_processed.rda"))




# Define some helper functions

# Load data from a file in a remote url or a local file appropriately
load_url <- function(url, cleanup = FALSE, ...) { 
  # Test whether the user gave a url
  if ( grepl('^http://.*$', url) ) { 
    file <- tempfile()
    download.file(url, destfile = file)
  } else { 
    file <- url
  }
  load(file, envir = parent.frame(), ...)
}

# Merge two summary data.frames
merge_branches <- function(upper, lower) { 
  data.frame(branch = c(rep('Upper branch', nrow(upper)), 
                        rep('Lower branch', nrow(lower))),
             rbind(upper, lower))
}

# Get the closest value in a vector
closest_to <- function(val, X, quiet = FALSE) { 
  index <- which( abs(val - X) == min(abs(val - X)) )
  new_value <- X[min(index)]
  if (!quiet) { 
    cat('Picked value ', new_value, ' (error: ', new_value - val,')\n', sep ='')
  }
  return(new_value)
}


no_legend <- theme(legend.position = "none")

# Function that produces a ggplot object representing the cross section
#   It requires ggplot2 > 2.0
do_cs_plot <- function(data, x, y, ytitle, main_title) { 
  
  ggplot(data) + 
    geom_line(aes_(x = substitute(x), y = substitute(y), 
                   group = ~branch, linetype = ~branch)) + 
    xlab(substitute(x)) + 
    ylab(ytitle) + 
    ggtitle(main_title) + 
    theme(plot.title = element_text(size = 14, hjust = 0, vjust = 0.5), 
          legend.position = c(1,1),
          legend.direction = 'horizontal', 
          legend.justification = c(1, 1), 
          legend.title = element_blank())
}

# Merge the plots in a three panes layout
merge_plots <- function(phasediag, cs_high, cs_mid, cs_low) { 
  grid.arrange(phasediag, 
               cs_high, cs_mid, cs_low, 
               ncol = 1, 
               heights = c(.5, .5/3, .5/3, .5/3))
}

# Values common for all models 
setwd(working_directory)

# Graphic output 
GRAPH_WIDTH <- 5
GRAPH_HEIGHT <- GRAPH_WIDTH * 2 # The 2x value ensures that 2D planes are squareish
FILE_PREFIX <- 'bifurc_diagram_'





# Musselbed model -------------------------

# Load data. We are interested in the summary data only. We discard the rest 
#   as RAM is limited.
load_url(files[['musselbed']], verbose = TRUE)
dat.mussel <- merge_branches(upper_branch[['DatBif']], 
                             lower_branch[['DatBif']])
dat.mussel <- subset(dat.mussel, delta > 0)
rm(upper_branch, lower_branch)

# Values for individual cross_sections
cross_section_high <- closest_to(0.60, dat.mussel[ ,'d'])
cross_section_mid  <- closest_to(0.24, dat.mussel[ ,'d'])
cross_section_low  <- closest_to(0, dat.mussel[ ,'d'])

# Data-specific values we add on the graph by hand
# NB: 0 and 1 contours are not plotted so no label is added (NA in values)
contour_breaks <- data.frame(breaks = c(0,    0.1,   0.2,  0.5,    0.8,  1),
                             lbl.x  = c(NA,   .135,  .125, .055, .024, NA),
                             lbl.y  = c(NA,   .30,  .256,  .20,   .09, NA),
                          lbl.angle = c(NA,   -45,   -45,   -45,    -55,  NA))

mussel_plane_plot <- 
  ggplot(NULL, aes(x = delta,y = d)) + 
    geom_raster(data = subset(dat.mussel, branch == "Lower branch" & 
                                          mean_cover_. == 0),
                fill = 'grey80') + 
    geom_contour(aes(z = mean_cover_., color = mean_cover_.), 
                 data = subset(dat.mussel, branch == "Upper branch"), 
                 color = 'black', breaks = contour_breaks[ ,'breaks']) + 
    geom_text(aes(x = lbl.x, 
                  y = lbl.y, 
                  angle = lbl.angle,
                  label = as.character(breaks)),
               data = contour_breaks) + 
    annotate("text", x = 0.0, y = cross_section_high, label = "a ▸", 
             hjust = 1, vjust = .5) + 
    annotate("text", x = 0.0, y = cross_section_mid, label = "b ▸", 
             hjust = 1, vjust = .5) + 
    annotate("text", x = 0.0, y = cross_section_low, label = "c ▸", 
             hjust = 1, vjust = .5) + 
    theme_minimal() + 
    ggtitle('Musselbed model')

# Build cross sections plots
mussel_high_cs_plot <- do_cs_plot(subset(dat.mussel, d == cross_section_high), 
                                  x = delta, y = mean_cover_., 
                                  ytitle = 'Mussel density', 
                                  main_title = 'a.')

mussel_mid_cs_plot <- do_cs_plot(subset(dat.mussel, d == cross_section_mid), 
                                 x = delta, y = mean_cover_., 
                                 ytitle = 'Mussel density', 
                                 main_title = 'b.') + 
                                 no_legend

mussel_low_cs_plot <- do_cs_plot(subset(dat.mussel, d == cross_section_low), 
                                 x = delta, y = mean_cover_., 
                                 ytitle = 'Mussel density', 
                                 main_title = 'c.') + 
                                 no_legend

# Merge plots and save
mussel_merged_plots <- merge_plots(mussel_plane_plot, 
                                   mussel_high_cs_plot,
                                   mussel_mid_cs_plot,
                                   mussel_low_cs_plot)

ggsave(mussel_merged_plots, 
       filename = paste0(output_figure_path, FILE_PREFIX, "musselbed.png"), 
       width = GRAPH_WIDTH, height = GRAPH_HEIGHT)




# Grazing model -------------------------

# Load data. We are interested in the summary data only. We discard the rest 
#   as RAM is limited.
load_url(files[['grazing']], verbose = TRUE)
dat.grazing <- merge_branches(upper_branch[['DatBif']], 
                             lower_branch[['DatBif']])
rm(upper_branch, lower_branch)

# Values for individual cross_sections
cross_section_high <- closest_to(0.18, dat.grazing[ ,'g'])
cross_section_mid  <- closest_to(0.072, dat.grazing[ ,'g'])
cross_section_low  <- closest_to(0, dat.grazing[ ,'g'])

# Data-specific values we add on the graph by hand
# NB: 0 and 1 contours are not plotted so no label is added (NA in values)
contour_breaks <- data.frame(breaks = c(0,    0.1,   0.5,    0.8,   1),
                             lbl.x  = c(NA,   .16,  .095,     0.030,   NA),
                             lbl.y  = c(NA,   .05,  .10,     0.1,  NA),
                         lbl.angle =  c(NA,   -45,    -55,   -73,  NA))

grazing_plane_plot <- 
  ggplot(NULL, aes(x = m0,y = g)) + 
    geom_raster(data = subset(dat.grazing, branch == "Lower branch" & 
                                          mean_cover_. == 0),
                fill = 'grey80') + 
    geom_contour(aes(z = mean_cover_., color = mean_cover_.), 
                 data = subset(dat.grazing, branch == "Upper branch"), 
                 color = 'black', breaks = contour_breaks[ ,'breaks']) + 
    geom_text(aes(x = lbl.x, 
                  y = lbl.y, 
                  angle = lbl.angle,
                  label = as.character(breaks)),
               data = contour_breaks) + 
    annotate("text", x = 0.0, y = cross_section_high, label = "a ▸", 
             hjust = 1, vjust = .5) + 
    annotate("text", x = 0.0, y = cross_section_mid, label = "b ▸", 
             hjust = 1, vjust = .5) + 
    annotate("text", x = 0.0, y = cross_section_low, label = "c ▸", 
             hjust = 1, vjust = .5) + 
    theme_minimal() + 
    ggtitle('Grazing model')

# Build cross sections plots
grazing_high_cs_plot <- do_cs_plot(subset(dat.grazing, g == cross_section_high), 
                                   x = m0, y = mean_cover_., 
                                   ytitle = 'Forest density', 
                                   main_title = 'a.')

grazing_mid_cs_plot <- do_cs_plot(subset(dat.grazing, g == cross_section_mid), 
                                  x = m0, y = mean_cover_., 
                                  ytitle = 'Forest density', 
                                  main_title = 'b.') + 
                                  no_legend
                                  
grazing_low_cs_plot <- do_cs_plot(subset(dat.grazing, g == cross_section_low), 
                                  x = m0, y = mean_cover_., 
                                  ytitle = 'Forest density', 
                                  main_title = 'c.') +
                                  no_legend


grazing_merged_plots <- merge_plots(grazing_plane_plot, 
                                    grazing_high_cs_plot,
                                    grazing_mid_cs_plot,
                                    grazing_low_cs_plot)


ggsave(grazing_merged_plots, 
       filename = paste0(output_figure_path, FILE_PREFIX, "grazing.png"), 
       width = GRAPH_WIDTH, height = GRAPH_HEIGHT)




# Forest Gap model -------------------------

# Load data. We are interested in the summary data only. We discard the rest 
#   as RAM is limited.
load_url(files[['forestgap']], verbose = TRUE)
dat.forestgap <- merge_branches(upper_branch[['DatBif']], 
                                lower_branch[['DatBif']])
rm(upper_branch, lower_branch)

# Values for individual cross_sections
cross_section_high <- closest_to(0.189, dat.forestgap[ ,'delta'])
cross_section_mid  <- closest_to(0.07, dat.forestgap[ ,'delta'])
cross_section_low  <- closest_to(0, dat.forestgap[ ,'delta'])

# Data-specific values we add on the graph by hand
# NB: 0 and 1 contours are not plotted so no label is added (NA in values)
contour_breaks <- data.frame(breaks = c(0,    0.1,     0.5,    0.8,   1),
                             lbl.x  = c(NA,   0.128,  0.048,  0.018,   NA),
                             lbl.y  = c(NA,   .05,     .1,     0.1,  NA),
                         lbl.angle =  c(NA,   -45,     -60,   -75,  NA))

forestgap_plane_plot <- 
  ggplot(NULL, aes(x = d, y = delta)) + 
    geom_raster(data = subset(dat.forestgap, branch == "Lower branch" & 
                                             mean_cover_. == 0),
                fill = 'grey80') + 
    geom_contour(aes(z = mean_cover_., color = mean_cover_.), 
                 data = subset(dat.forestgap, branch == "Upper branch"), 
                 color = 'black', breaks = contour_breaks[ ,'breaks']) + 
    geom_text(aes(x = lbl.x, 
                  y = lbl.y, 
                  angle = lbl.angle,
                  label = as.character(breaks)),
               data = contour_breaks) + 
    annotate("text", x = 0.0, y = cross_section_high, label = "a ▸", 
             hjust = 1, vjust = .5) + 
    annotate("text", x = 0.0, y = cross_section_mid, label = "b ▸", 
             hjust = 1, vjust = .5) + 
    annotate("text", x = 0.0, y = cross_section_low, label = "c ▸", 
             hjust = 1, vjust = .5) + 
    theme_minimal() + 
    ggtitle('Forest Gap model')

# Build cross sections plots
forestgap_high_cs_plot <- do_cs_plot(subset(dat.forestgap, delta == cross_section_high), 
                                   x = d, y = mean_cover_., 
                                   ytitle = 'Forest density', 
                                   main_title = 'a.')

forestgap_mid_cs_plot <- do_cs_plot(subset(dat.forestgap, delta == cross_section_mid), 
                                  x = d, y = mean_cover_., 
                                  ytitle = 'Forest density', 
                                  main_title = 'b.') + 
                                  no_legend
                                  
forestgap_low_cs_plot <- do_cs_plot(subset(dat.forestgap, delta == cross_section_low), 
                                  x = d, y = mean_cover_., 
                                  ytitle = 'Forest density', 
                                  main_title = 'c.') +
                                  no_legend

forestgap_merged_plots <- merge_plots(forestgap_plane_plot, 
                                      forestgap_high_cs_plot,
                                      forestgap_mid_cs_plot,
                                      forestgap_low_cs_plot)

ggsave(forestgap_merged_plots, 
       filename = paste0(output_figure_path, FILE_PREFIX, "forestgap.png"), 
       width = GRAPH_WIDTH, height = GRAPH_HEIGHT)






# All plots on the same plate --------------------------

all_plots_merged <- grid.arrange(mussel_plane_plot, 
                                 grazing_plane_plot, 
                                 forestgap_plane_plot,
                                 
                                 mussel_high_cs_plot, 
                                 grazing_high_cs_plot, 
                                 forestgap_high_cs_plot, 
                                 
                                 mussel_mid_cs_plot, 
                                 grazing_mid_cs_plot, 
                                 forestgap_mid_cs_plot, 
                                 
                                 mussel_low_cs_plot, 
                                 grazing_low_cs_plot, 
                                 forestgap_low_cs_plot, 
                                 ncol = 3,
                                 heights = c( .5, .5/3, .5/3, .5/3),
                                 widths  = c(1/3, 1/3, 1/3))

ggsave(all_plots_merged, 
       filename = paste0(output_figure_path, FILE_PREFIX, "all.png"), 
       width = GRAPH_WIDTH * 3, 
       height = GRAPH_HEIGHT)

ggsave(all_plots_merged, 
       filename = paste0(output_figure_path, FILE_PREFIX, "all.pdf"), 
       width = GRAPH_WIDTH * 3, 
       height = GRAPH_HEIGHT)


# Done.
