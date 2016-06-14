# 
# 
# This file contains code that computes the temporal indicators
# 

# State dependencies
library(dplyr)
library(plyr)
library(tidyr)
library(moments)
library(ggplot2)

# Define output path 
# /!\ NB: folders must have a trailing /
working_directory  <- '/home/alex/work/2014-2015/SpatialStress/ifcam_output/figures/'
output_figure_path <- "/home/alex/work/2015-2016/SpatialStress/ifcam_output/figures/temporal_indicators/" # same folder

# Define data paths 
# /!\ NB: folders must have a trailing /
data_folder <- '/home/alex/work/2014-2015/SpatialStress/ifcam_output/data/'
# Use this folder if you do not have the results on your computer: it will fetch
#   the latest ones.
# data_folder <- "http://alex.lecairn.org/ifcam/" 
files <- list(musselbed = paste0(data_folder, "result_musselbed_cs_processed.rda"),
              grazing   = paste0(data_folder, "result_grazing_cs_processed.rda"),
              forestgap = paste0(data_folder, "result_forestgap_cs_processed.rda"))


# Graphic output 
GRAPH_WIDTH <- 9
GRAPH_HEIGHT <- GRAPH_WIDTH  # The 2x value ensures that 2D planes are squareish
FILE_PREFIX  <- paste0('temporal_indicators_')
output_figure_path <- "/home/alex/work/2014-2015/SpatialStress/ifcam_output/figures/temporal_indicators/" # same folder



# Computation parameters -------------------------------
# We will use the LASTN values to compute indicators
LASTN <- 300



# Helper functions -------------------------

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

# Merge two summary data.frames with a key describing the origin of the data
merge_df <- function(df1, df2, key_origin1, key_origin2, key_colname) { 
  # Create key column
  key_col <- data.frame(c(rep(key_origin1, nrow(df1)), 
                          rep(key_origin2, nrow(df2))), 
                        stringsAsFactors = FALSE)
  
  new_df <- data.frame(key_col, rbind(df1, df2))
  names(new_df) <- c(key_colname, names(df1))
  
  return(new_df)
}

# Repair df names when it is set to X. X0 or X..1
repair_names <- function(df) { 
  old <- c('X\\.$', 'X0$', 'X\\.\\.1$')
  new <- c('plus', 'zero', 'minus')
  for ( i in seq_along(old) ) { 
    names(df)[grepl(old[i], names(df))] <- new[i]
  }
  return(df)
}

# Reformat time-series
merge_global_local <- function(list_of_both) { 
  # Add time column from 1 to nrow(df)
  list_of_both <- lapply(list_of_both, 
                         function(df) { data.frame(time = seq.int(nrow(df)), df) })
  merge_df(list_of_both[["global"]], list_of_both[["local"]], 
           key_origin1 = 'global', key_origin2 = "local", 
           key_colname = "density_type")
}

# Merge the time series into a data.frame  
merge_ts_into_df <- function(branch) { 
  # Massage data so everything fits into a data.frame
  ts <- lapply(branch[['time_series']], merge_global_local)
  ts <- Map(function(df, n) data.frame(ID = n, df), 
            ts, branch[['DatBif']][ ,'ID'])
  ts <- rbind_all(ts) # use dplyr's fast rbind instead of do.call(rbind, ts)
  ts <- repair_names(ts)
  ts
} 



# Temporal indicator functions -------------------------

ic_autocor_lag1 <- function(X) acf(X, plot = FALSE, lag.max = 1)[['acf']][2]
ic_variance     <- function(X) var(X) # X is a vector
ic_skewness     <- function(X) skewness(X) # X is a vector

ic_all <- function(df, col, replace_NA_by = NA) { 
  values <- ifelse(is.na(df[ ,col]), replace_NA_by, df[ ,col])
  rbind(data.frame(Indicator = 'Mean',      value = mean(values)), 
        data.frame(Indicator = 'lag-1 Autocorr.', value = ic_autocor_lag1(values)), 
        data.frame(Indicator = 'Skewness',  value = ic_skewness(values)),
        data.frame(Indicator = 'Variance',  value = ic_variance(values))) 
}




# Musselbed model -------------------------

# Reformat times-series into a big data frame
load(files[["musselbed"]], verbose = TRUE)

tseries.mussel <- merge_ts_into_df(branch)

# Get LASTN values
tseries.mussel <- ddply(tseries.mussel, ~ density_type + ID, 
                        tail, n = LASTN, 
                        .progress = 'time')

# Check if equilibrium is attained
# ggplot( subset(tseries.mussel, density_type == 'global') ) + 
#   geom_line(aes(time, plus, group = ID), alpha = .2)

# Compute indicators
ic.mussel <- ddply(tseries.mussel, ~ density_type + ID, 
                   ic_all, col = "plus", replace_NA_by = 0,
                   .progress = "time")

# Join with data.frame containing parameters
ic.mussel <- join(branch[["DatBif"]], ic.mussel, by = 'ID')

# Create plot
plot.mussel <- 
  ggplot( subset(ic.mussel, density_type == "global") ) + 
    geom_line(aes(x = delta, y = value)) + 
    facet_grid(Indicator ~ d, scales = "free_y", labeller = label_both)

ggsave(paste0(output_figure_path, FILE_PREFIX, "musselbed.pdf"),
       width = GRAPH_WIDTH, height = GRAPH_HEIGHT,
       plot = plot.mussel)




# Grazing model -------------------------

# Reformat times-series into a big data frame
load(files[["grazing"]], verbose = TRUE)

tseries.grazing <- merge_ts_into_df(branch)

# Get LASTN values
tseries.grazing <- ddply(tseries.grazing, ~ density_type + ID, 
                         tail, n = LASTN, 
                         .progress = 'time')

# Check if equilibrium is attained
ggplot( subset(tseries.grazing, density_type == 'global') ) + 
  geom_line(aes(time, plus, group = ID), alpha = .2)


# Compute indicators
ic.grazing <- ddply(tseries.grazing, ~ density_type + ID, 
                   ic_all, col = "plus", replace_NA_by = 0,
                   .progress = "time")

# Join with data.frame containing parameters
ic.grazing <- join(branch[["DatBif"]], ic.grazing, by = 'ID')

# Create plot
plot.grazing <- 
  ggplot( subset(ic.grazing, density_type == "global") ) + 
    geom_line(aes(x = m0, y = value)) + 
    facet_grid(Indicator ~ g, scales = "free_y", labeller = label_both)

ggsave(paste0(output_figure_path, FILE_PREFIX, "grazing.pdf"),
       width = GRAPH_WIDTH, height = GRAPH_HEIGHT,
       plot = plot.grazing)




# Forestgap model -------------------------

# Reformat times-series into a big data frame
load(files[["forestgap"]], verbose = TRUE)

tseries.forestgap <- merge_ts_into_df(branch)

# Get LASTN values
tseries.forestgap <- ddply(tseries.forestgap, ~ density_type + ID, 
                           tail, n = LASTN, 
                           .progress = 'time')

# Check if equilibrium is attained
ggplot( subset(tseries.forestgap, density_type == 'global') ) + 
  geom_line(aes(time, plus, group = ID), alpha = .2)

# Compute indicators
ic.forestgap <- ddply(tseries.forestgap, ~ density_type + ID, 
                      ic_all, col = "plus", replace_NA_by = 0,
                      .progress = "time")

# Join with data.frame containing parameters
ic.forestgap <- join(branch[["DatBif"]], ic.forestgap, by = 'ID')

# Create plot
plot.forestgap <- 
  ggplot( subset(ic.forestgap, density_type == "global") ) + 
    geom_line(aes(x = d, y = value)) + 
    facet_grid(Indicator ~ delta, scales = "free_y", labeller = label_both)

ggsave(paste0(output_figure_path, FILE_PREFIX, "forestgap.pdf"),
       width = GRAPH_WIDTH, height = GRAPH_HEIGHT,
       plot = plot.forestgap)

