# 
# 
# This file contains code that plots spectral indicators along the stress 
#   gradient. Spectral indicators include the r-spectrum and SDR that 
#   derives from it.
# 

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
GRAPH_HEIGHT <- GRAPH_WIDTH  # The 2x value ensures that 2D planes are squareish
FILE_PREFIX  <- paste0('spectral_indicators_',CG_SUBSIZE,"x",CG_SUBSIZE,"_")
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


load_url(files[['grazing_cs']])

id_subset <- with(branch, runif(nrow(DatBif)) < .05) # 5%

result <- add_spectral_indic(branch[['DatBif']][id_subset, ], 
                             branch[['snaps']][id_subset], 
                             NULL_REPS, SDR_LOW_RANGE, SDR_HIGH_RANGE)





