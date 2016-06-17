# 
# This file produces a figure with generic indicator trends along cross 
#   sections. 
# 
# This code makes use of the package functions. 
# 

library(ggplot2)
library(plyr) 
library(grid) 
library(gridExtra) # for grid.arrange
library(doParallel) # for grid.arrange

# Install latest spatialwarnings
library(spatialwarnings) 
# install_github('fdschneider/spatial_warnings') # mind the _

# Indicator computation parms
CG_SUBSIZE <- 10
NULL_REPS <- 99

# Define output path 
# /!\ NB: folders must have a trailing /
working_directory  <- '/home/alex/work/2014-2015/SpatialStress/ifcam_output/figures/'
setwd(working_directory)


# Graphic output 
GRAPH_WIDTH <- 9
GRAPH_HEIGHT <- GRAPH_WIDTH  # The 2x value ensures that 2D planes are squareish
FILE_PREFIX  <- paste0('generic_indicators_',CG_SUBSIZE,"x",CG_SUBSIZE,"_")
output_figure_path <- "/home/alex/work/2014-2015/SpatialStress/ifcam_output/figures/generic_indicators/" # same folder

# Define data paths 
# /!\ NB: folders must have a trailing /
data_folder <- '/home/alex/work/2014-2015/SpatialStress/ifcam_output/data/'
# Use this folder if you do not have the results on your computer: it will fetch
#   the latest ones from Alex's server. Do not forget to use the proxy if 
#   needed (Sys.setenv(http_proxy=...)).
# data_folder <- "http://alex.lecairn.org/ifcam/"
files <- list(musselbed = paste0(data_folder, "result_musselbed_processed.rda"),
              grazing   = paste0(data_folder, "result_grazing_processed.rda"),
              forestgap = paste0(data_folder, "result_forestgap_processed.rda"),
              grazing_cs = paste0(data_folder, "result_grazing_cs_processed.rda"),
              forestgap_cs = paste0(data_folder, "result_forestgap_cs_processed.rda"),
              musselbed_cs = paste0(data_folder, "result_musselbed_cs_processed.rda"))

# Use parallelism ? 
PARALLEL <- FALSE

# Register parallel backend
if (PARALLEL) { 
  registerDoParallel(cores = 10)
} else { 
  registerDoSEQ()
}



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
  gc()
  return(TRUE)
}

add_generic_indic <- function(df, # data.frame of all simulations (DatBif)
                              matrices, # corresponding matrices (list of lists)
                              subsize, # coarse-graining subsize
                              null_replicates, # number of null replicates
                              ...) { 
  
  # Compute generic indicators (uses summary method for generic spews)
  result <- llply(matrices, generic_spews, subsize = subsize, ...)
  result <- llply(result, indictest, null_replicates = null_replicates, ...)
  
  # Add ID column for merge
  ids <- Map(rep, df[ ,'ID'], lapply(result, nrow))
  ids <- unlist(ids)
  result <- data.frame(ID = ids, do.call(rbind, result))
  
  # Merge all results in one df and return
  merged <- join(df, result, by = 'ID')
  
  # Clean up mem as soon as possible
  rm(df, matrices)
  gc()
  
  return(merged)
}




# Grazing model
# -------------------------

load_url(files[['grazing_cs']], verbose = TRUE)

upper_result <- add_generic_indic(branch[['DatBif']], 
                                  branch[['snaps']], 
                                  subsize = CG_SUBSIZE, 
                                  null_replicates = NULL_REPS, 
                                  .progress = 'time',
                                  .parallel = PARALLEL)

upper_result_trend <- ddply(upper_result, ~ ID + indicator + g + m0, 
                            summarise, 
                            value.mean = mean(value), 
                            value.sd_above = mean(value) + sd(value), 
                            value.sd_below = mean(value) - sd(value),
                            null.mean = mean(null_mean), 
                            null.sd_above = mean(null_mean) + sd(null_mean), 
                            null.sd_below = mean(null_mean) - sd(null_mean),
                            .progress = 'time')

# Plot chosen cross-sections
grazing_csplot <- 
  ggplot(upper_result_trend) + 
    geom_ribbon(aes(m0, ymin = null.sd_below, ymax = null.sd_above), alpha = .2, 
                color = 'red') + 
    geom_line(aes(m0, null.mean), color = 'red') + 
    geom_ribbon(aes(m0, ymin = value.sd_below, ymax = value.sd_above), 
                alpha = .4, color = 'black') + 
    geom_line(aes(m0, value.mean), alpha = .6, color = 'black')  + 
    facet_grid( indicator ~ g, scales = 'free_y', labeller = label_both) + 
    ggtitle('Grazing model')

ggsave(paste0(output_figure_path, FILE_PREFIX, "grazing.pdf"),
       width = GRAPH_WIDTH, height = GRAPH_HEIGHT,
       plot = grazing_csplot)

rm(branch)
gc()

# # Explore some images
# jj <- subset(upper_result_trend, m0 == closest_to(0.008, m0) & 
#                                    g %in% c(0, 
#                                             closest_to(0.072, g), 
#                                             closest_to(0.18, g)))
# ids <- unique(jj$ID)
# 
# for (id in ids) { 
#   dev.new()
#   print(lattice::levelplot(upper_branch[["snaps"]][[id]][[1]]))
# }





# Forestgap model
# -------------------------

load_url(files[['forestgap_cs']], verbose = TRUE)

upper_result <- add_generic_indic(branch[['DatBif']], 
                                  branch[['snaps']], 
                                  subsize = CG_SUBSIZE, 
                                  null_replicates = NULL_REPS, 
                                  .progress = 'time',
                                  .parallel = PARALLEL)

upper_result_trend <- ddply(upper_result, ~ ID + indicator + d + delta, 
                            summarise, 
                            value.mean = mean(value), 
                            value.sd_above = mean(value) + sd(value), 
                            value.sd_below = mean(value) - sd(value),
                            null.mean = mean(null_mean), 
                            null.sd_above = mean(null_mean) + sd(null_mean), 
                            null.sd_below = mean(null_mean) - sd(null_mean),
                            .progress = 'time')

forestgap_csplot <- 
  ggplot(upper_result_trend) + 
    geom_ribbon(aes(d, ymin = null.sd_below, ymax = null.sd_above), alpha = .2, 
                color = 'red') + 
    geom_line(aes(d, null.mean), color = 'red') + 
    geom_ribbon(aes(d, ymin = value.sd_below, ymax = value.sd_above), 
                alpha = .4, color = 'black') + 
    geom_line(aes(d, value.mean), alpha = .6, color = 'black')  + 
    facet_grid( indicator ~ delta, scales = 'free_y', labeller = label_both) + 
    ggtitle('Forestgap model')

ggsave(paste0(output_figure_path, FILE_PREFIX, "forestgap.pdf"),
       width = GRAPH_WIDTH, height = GRAPH_HEIGHT,
       plot = forestgap_csplot)

rm(branch)
gc()




# Musselbed model
# -------------------------

load_url(files[['musselbed_cs']], verbose = TRUE)

upper_result <- add_generic_indic(branch[['DatBif']], 
                                  branch[['snaps']], 
                                  subsize = CG_SUBSIZE, 
                                  null_replicates = NULL_REPS, 
                                  .progress = 'time',
                                  .parallel = PARALLEL)

upper_result_trend <- ddply(upper_result, ~ ID + indicator + d + delta, 
                            summarise, 
                            value.mean = mean(value), 
                            value.sd_above = mean(value) + sd(value), 
                            value.sd_below = mean(value) - sd(value),
                            null.mean = mean(null_mean), 
                            null.sd_above = mean(null_mean) + sd(null_mean), 
                            null.sd_below = mean(null_mean) - sd(null_mean),
                            .progress = 'time')

musselbed_csplot <- 
  ggplot(upper_result_trend) + 
    geom_ribbon(aes(delta, ymin = null.sd_below, ymax = null.sd_above), alpha = .2, 
                color = 'red') + 
    geom_line(aes(delta, null.mean), color = 'red') + 
    geom_ribbon(aes(delta, ymin = value.sd_below, ymax = value.sd_above), 
                alpha = .4, color = 'black') + 
    geom_line(aes(delta, value.mean), alpha = .6, color = 'black')  + 
    facet_grid( indicator ~ d, scales = 'free_y', labeller = label_both) + 
    ggtitle('Musselbed model')

ggsave(paste0(output_figure_path, FILE_PREFIX, "musselbed.pdf"),
       width = GRAPH_WIDTH, height = GRAPH_HEIGHT,
       plot = musselbed_csplot)

rm(branch)
gc()



