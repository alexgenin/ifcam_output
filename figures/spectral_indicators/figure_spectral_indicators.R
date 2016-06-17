# 
# 
# This file contains code that plots spectral indicators along the stress 
#   gradient. Spectral indicators include the r-spectrum and SDR that 
#   derives from it.
# 

library(ggplot2)
library(plyr) 
library(grid) 
library(gridExtra) # for grid.arrange
library(doParallel) 
library(devtools) 

# Install latest spatialwarnings
install_github('fdschneider/spatial_warnings') # mind the _
library(spatialwarnings) 

# Indicator computation parms
NULL_REPS <- 4

# SDR-related variables
SDR_LOW_RANGE  <- c(0, .2) # lower 20%
SDR_HIGH_RANGE <- c(.8, 1) # higher 20%

# Define output path 
# /!\ NB: folders must have a trailing /
working_directory  <- '/home/alex/work/2014-2015/SpatialStress/ifcam_output/figures/'
setwd(working_directory)


# Graphic output 
GRAPH_WIDTH <- 9
GRAPH_HEIGHT <- 3  # The 2x value ensures that 2D planes are squareish
FILE_PREFIX  <- paste0('spectral_indicators_')
output_figure_path <- paste0(working_directory, 'spectral_indicators/')

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
PARALLEL <- TRUE

# Register parallel backend
if (PARALLEL) { 
  registerDoParallel(cores = 10)
} else { 
  registerDoSEQ()
}

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


# Function that adds the spectral indicator column to the DatBif df
add_spectral_indic <- function(df, # data.frame of all simulations (DatBif)
                               matrices, # corresponding matrices (list of lists)
                               null_replicates, # number of null replicates
                               low_range, high_range, 
                               ...) { 
  
  # Compute generic indicators (uses summary method for generic spews)
  result <- llply(matrices, spectral_spews,
                  sdr_low_range = low_range, 
                  sdr_high_range = high_range, 
                  .progress = 'time', 
                  .parallel = PARALLEL,  
                  ...)
  
  result <- llply(result, indictest, 
                  null_replicates = null_replicates, 
                  .progress = 'time',
                  .parallel = PARALLEL,
                  ...)
  
  # Add ID column for merge
  ids <- Map(rep, df[ ,'ID'], lapply(result, nrow))
  
  result <- data.frame(ID = unlist(ids), 
                       do.call(rbind, result))
  
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

result <- add_spectral_indic(branch[['DatBif']],
                             branch[['snaps']],
                             NULL_REPS, SDR_LOW_RANGE, SDR_HIGH_RANGE)

result_summ <- ddply(subset(result, type == 'sdr'), 
                     ~ ID + g + m0, 
                     summarise, 
                     value.mean = mean(value), 
                     value.sd_above = mean(value) + sd(value), 
                     value.sd_below = mean(value) - sd(value),
                     null.mean = mean(null_mean), 
                     null.sd_above = mean(null_mean) + sd(null_mean), 
                     null.sd_below = mean(null_mean) - sd(null_mean),
                     .progress = 'time')

grazing_sdr_plot <- 
  ggplot(result_summ) + 
    geom_ribbon(aes(m0, ymin = null.sd_below, ymax = null.sd_above), alpha = .6, 
                fill = 'red') + 
    geom_line(aes(m0, null.mean), color = 'red') + 
    geom_ribbon(aes(m0, ymin = value.sd_below, ymax = value.sd_above), 
                alpha = .4, fill = 'black') + 
    geom_line(aes(m0, value.mean), alpha = .6, color = 'black')  + 
    facet_grid( . ~ g, scales = 'free_y', labeller = label_both) + 
    ylab('Spectral Density Ratio') + 
    ggtitle('Grazing model')

ggsave(paste0(output_figure_path, FILE_PREFIX, "grazing.pdf"),
       width = GRAPH_WIDTH, height = GRAPH_HEIGHT,
       plot = grazing_sdr_plot)



# Musselbed model
# -------------------------

load_url(files[['musselbed_cs']], verbose = TRUE)

result <- add_spectral_indic(branch[['DatBif']],
                             branch[['snaps']],
                             NULL_REPS, SDR_LOW_RANGE, SDR_HIGH_RANGE)


result_summ <- ddply(subset(result, type == 'sdr'), 
                     ~ ID + d + delta, 
                     summarise, 
                     value.mean = mean(value), 
                     value.sd_above = mean(value) + sd(value), 
                     value.sd_below = mean(value) - sd(value),
                     null.mean = mean(null_mean), 
                     null.sd_above = mean(null_mean) + sd(null_mean), 
                     null.sd_below = mean(null_mean) - sd(null_mean),
                     .progress = 'time')

musselbed_sdr_plot <- 
  ggplot(result_summ) + 
    geom_ribbon(aes(delta, ymin = null.sd_below, ymax = null.sd_above), alpha = .6, 
                fill = 'red') + 
    geom_line(aes(delta, null.mean), color = 'red') + 
    geom_ribbon(aes(delta, ymin = value.sd_below, ymax = value.sd_above), 
                alpha = .4, fill = 'black') + 
    geom_line(aes(delta, value.mean), alpha = .6, color = 'black')  + 
    facet_grid( . ~ d, scales = 'free_y', labeller = label_both) + 
    ylab('Spectral Density Ratio') + 
    ggtitle('Musselbed model')

ggsave(paste0(output_figure_path, FILE_PREFIX, "musselbed.pdf"),
       width = GRAPH_WIDTH, height = GRAPH_HEIGHT,
       plot = musselbed_sdr_plot)



# Forestgap model
# -------------------------

load_url(files[['forestgap_cs']], verbose = TRUE)

result <- add_spectral_indic(branch[['DatBif']],
                             branch[['snaps']],
                             NULL_REPS, SDR_LOW_RANGE, SDR_HIGH_RANGE)

result_summ <- ddply(subset(result, type == 'sdr'), 
                     ~ ID + d + delta, 
                     summarise, 
                     value.mean = mean(value), 
                     value.sd_above = mean(value) + sd(value), 
                     value.sd_below = mean(value) - sd(value),
                     null.mean = mean(null_mean), 
                     null.sd_above = mean(null_mean) + sd(null_mean), 
                     null.sd_below = mean(null_mean) - sd(null_mean),
                     .progress = 'time')

forestgap_sdr_plot <- 
  ggplot(result_summ) + 
    geom_ribbon(aes(d, ymin = null.sd_below, ymax = null.sd_above), alpha = .6, 
                fill = 'red') + 
    geom_line(aes(d, null.mean), color = 'red') + 
    geom_ribbon(aes(d, ymin = value.sd_below, ymax = value.sd_above), 
                alpha = .4, fill = 'black') + 
    geom_line(aes(d, value.mean), alpha = .6, color = 'black')  + 
    facet_grid( . ~ delta, scales = 'free_y', labeller = label_both) + 
    ylab('Spectral Density Ratio') + 
    ggtitle('Forestgap model')

ggsave(paste0(output_figure_path, FILE_PREFIX, "forestgap.pdf"),
       width = GRAPH_WIDTH, height = GRAPH_HEIGHT,
       plot = forestgap_sdr_plot)


