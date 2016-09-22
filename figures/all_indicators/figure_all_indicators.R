# 
# 
# This file plots for each model a big figure with all indicators, with three 
#   columns (one for each value of heterogeneous stressor) and 8 to 10 columns
# 


library(ggplot2)
library(ggthemes)
library(plyr) 
library(grid) 
library(gridExtra) # for grid.arrange
library(doParallel) 
library(devtools) 
library(moments) 

# Install latest spatialwarnings
install_github('fdschneider/spatial_warnings') # mind the _
library(spatialwarnings) 

# Indicator computation parms
CG_SUBSIZE <- 10
LASTN <- 300 # for temporal indicators, take the last xx time steps
NULL_REPS <- 499

# Define output path 
# /!\ NB: folders must have a trailing /
working_directory  <- '/home/alex/work/2014-2015/SpatialStress/ifcam_output/figures/'
setwd(working_directory)

# Graphic output 
GRAPH_WIDTH <- 8
GRAPH_HEIGHT <- GRAPH_WIDTH * (3/2)
FILE_PREFIX  <- paste0('all_indicators_',CG_SUBSIZE,"x",CG_SUBSIZE,"_")
output_figure_path <- "./all_indicators/" 

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
  # Note: we explicitely request a PSOCK cluster as fork()-based method crash...
  cl <- makePSOCKcluster(3)
  registerDoParallel(cl = cl)
} else { 
  registerDoSEQ()
}


# Graph theme
theme_ifcam <- 
  theme_bw() + 
  theme(axis.title.y = element_blank(), 
        strip.background = element_blank(), 
        panel.border = element_rect(fill = NA, linetype = 'dotted', color = 'grey20'))


# Define some helper functions

# Temporal indicator functions -------------------------

ic_autocor_lag1 <- function(X) acf(X, plot = FALSE, lag.max = 1)[['acf']][2]
ic_variance     <- function(X) var(X) # X is a vector
ic_skewness     <- function(X) skewness(X) # X is a vector
ic_temporal_sdr <- function(X) { 
  spec <- spectrum(X, plot = FALSE)[['spec']]
  if ( any(spec > 0) ) { 
    sdr <- mean(spec[seq.int(1, round(length(spec)) * .2)]) / 
            mean(spec[seq.int(round(length(spec)) * .8, length(spec))])
  } else { 
    sdr <- NA
  }
  return(sdr)
}

temp_ic_all <- function(df, col, replace_NA_by = NA) { 
  values <- ifelse(is.na(df[ ,col]), replace_NA_by, df[ ,col])
  values <- tail(values, LASTN)
  rbind(data.frame(indicator = 'lag-1 Autocorr.',   value = ic_autocor_lag1(values)), 
        data.frame(indicator = 'Temporal Skewness', value = ic_skewness(values)),
        data.frame(indicator = 'Temporal Variance', value = ic_variance(values)), 
        data.frame(indicator = 'Temporal sdr',      value = ic_temporal_sdr(values)))
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

add_indics <- function(df, # data.frame of all simulations (DatBif)
                       matrices, # corresponding matrices (list of lists)
                       timeseries,
                       subsize, # coarse-graining subsize
                       null_replicates, # number of null replicates
                       ...) { 
  gc()
  
  # Compute spatial spectral indicators 
  result_specic <- llply(matrices, spectral_spews, 
                         sdr_low_range  = c(0, .2),      # lower 20%
                         sdr_high_range = c(.8, 1), ...) # upper 20%
  result_specic <- llply(result_specic, indictest, 
                         null_replicates = null_replicates,
                         ...)
  result_specic <- Map(function(n,df) { data.frame(ID = n, df) }, 
                       df[ ,"ID"], result_specic)
  result_specic <- do.call(rbind, result_specic)
  
  # Compute spatial generic indicators and null values
  result_genic <- llply(matrices, generic_spews, subsize = subsize, 
                        detrend = FALSE, moranI_coarse_grain = FALSE, ...)
  result_genic <- llply(result_genic, indictest, 
                        null_replicates = null_replicates, 
                        ...)
  result_genic <- Map(function(n,df) { data.frame(ID = n, df) }, 
                       df[ ,"ID"], result_genic)
  result_genic <- do.call(rbind, result_genic)
  
  
  # Merge spatial sdr and spatial genindic into one df
  # Drop spectrum-related data as it will blow up the memory
  result_specic <- result_specic[result_specic[ ,'type'] == "sdr", ]
  result_specic[ ,'indicator'] <- result_specic[ ,'type']
  cols <- c('ID', 'replicate', 'indicator', 'value', 'null_mean')
  result_spatial <- rbind(result_specic[ , cols], result_genic[ , cols])
  
  # Compute temporal indicators
  result_ic <- llply(timeseries, function(elem) temp_ic_all(elem[['global']], "+") )
  result_ic <- Map(function(n,e) { c(ID = n, e) }, df[ ,'ID'], result_ic)
  result_ic <- do.call(rbind, lapply(result_ic, as.data.frame))
  
  # Merge all results in one df and return
  result_all <- rbind.fill(result_spatial, result_ic)
  merged <- join(df, result_all, by = "ID")
  
  # Clean up mem as soon as possible
  rm(result_ic, result_specic, result_genic)
  gc()
  
  return(merged)
}

summarise_trend <- function(df, formula, ...) { 
  
  upper_result_trend <- ddply(df, formula, summarise, 
        value.mean = mean(value), 
        value.sd_above = mean(value) + sd(value), 
        value.sd_below = mean(value) - sd(value),
        null.mean = mean(null_mean), 
        null.sd_above = mean(null_mean) + sd(null_mean), 
        null.sd_below = mean(null_mean) - sd(null_mean),
        ...)
  
  # Massage data
  #   - Reorder and rename factors
  upper_result_trend[ ,'indicator'] <- 
    factor(upper_result_trend[ ,'indicator'], 
          levels = c('Mean', 'Moran\'s I', 'sdr', 'Variance', 'Skewness', 
                      'lag-1 Autocorr.', 'Temporal sdr', 'Temporal Variance', 
                      'Temporal Skewness'))
  levels(upper_result_trend[ ,'indicator']) <- 
    c('Mean', 'Moran\'s I', 'Spatial SDR', 'Spatial Variance', 'Spatial Skewness', 
      'lag-1 Autocorr.', 'Temporal SDR', 'Temporal Variance', 'Temporal Skewness')
  
  #   - Keep only data points with \rho_+ > 0
  upper_result_trend <- subset(upper_result_trend, 
                               ID %in% ID[indicator == "Mean" & value.mean > 0])
  #   - Drop mean 
  upper_result_trend <- subset(upper_result_trend, indicator != "Mean")
  
  upper_result_trend
}

# Find and give the tipping point coordinates 
create_shift_data_frame <- function(trends, heteroparm, homoparm) { 
  ddply(trends, heteroparm, 
        function(df) { data.frame(xshift = max(df[ ,homoparm])) })
}





# Grazing model
# -------------------------
if ( !exists("branch") ) { 
  load_url(files[['grazing_cs']], verbose = TRUE)
}

upper_result <- add_indics(branch[['DatBif']], 
                           branch[['snaps']], 
                           branch[['time_series']],
                           subsize = CG_SUBSIZE, 
                           null_replicates = NULL_REPS, 
                           .progress = 'time',
                           .parallel = PARALLEL)

upper_result_trend <- summarise_trend(upper_result, ~ ID + indicator + g + m0)

shift_df <- create_shift_data_frame(upper_result_trend, "g", "m0")

upper_result_trend[ ,'pretty_g'] <- with(upper_result_trend, paste0("g = ", g))
shift_df[ ,'pretty_g'] <- with(shift_df, paste0("g = ", g))

save(upper_result, upper_result_trend, shift_df, 
     file = paste0(output_figure_path, "grazing_all_indic_", 
                   CG_SUBSIZE, "x", CG_SUBSIZE, ".rda"))
     

# save.image('grazing_allindic.rda')

# Plot chosen cross-sections
grazing_csplot <- 
  ggplot(upper_result_trend) + 
    geom_ribbon(aes(m0, ymin = null.sd_below, ymax = null.sd_above), 
                alpha = .2, color = 'red') + 
    geom_line(aes(m0, null.mean), color = 'red') + 
    geom_ribbon(aes(m0, ymin = value.sd_below, ymax = value.sd_above), 
                alpha = .4, color = 'black') + 
    geom_line(aes(m0, value.mean), alpha = .9, color = 'black')  + 
    geom_vline(aes(xintercept = xshift), 
               color = 'grey50', linetype = 'dashed',
               data = shift_df) + 
    facet_grid(indicator ~ pretty_g, scales = 'free_y', 
               switch = "y") + 
    ggtitle('Grazing model') + 
    xlab( expression(m[0]) ) + 
    theme_ifcam

ggsave(paste0(output_figure_path, FILE_PREFIX, "grazing.pdf"),
       width = GRAPH_WIDTH, height = GRAPH_HEIGHT,
       plot = grazing_csplot)

rm(branch)
gc()





# Forestgap model
# -------------------------

if ( !exists("branch") ) { 
  load_url(files[['forestgap_cs']], verbose = TRUE)
}

upper_result <- add_indics(branch[['DatBif']], 
                           branch[['snaps']], 
                           branch[['time_series']],
                           subsize = CG_SUBSIZE, 
                           null_replicates = NULL_REPS, 
                           .progress = 'time',
                           .parallel = PARALLEL)
# save.image("forestgap.dump.rda")

upper_result_trend <- summarise_trend(upper_result, ~ ID + indicator + delta + d)

shift_df <- create_shift_data_frame(upper_result_trend, "delta", "d")

upper_result_trend[ ,'pretty_delta'] <- with(upper_result_trend, paste0("delta = ", delta))
shift_df[ ,'pretty_delta'] <- with(shift_df, paste0("delta = ", delta))

save(upper_result, upper_result_trend, shift_df, 
     file = paste0(output_figure_path, "forestgap_all_indic_", 
                   CG_SUBSIZE, "x", CG_SUBSIZE, ".rda"))
     

# save.image('grazing_allindic.rda')

# Plot chosen cross-sections
forestgap_csplot <- 
  ggplot(upper_result_trend) + 
    geom_ribbon(aes(d, ymin = null.sd_below, ymax = null.sd_above), 
                alpha = .2, color = 'red') + 
    geom_line(aes(d, null.mean), color = 'red') + 
    geom_ribbon(aes(d, ymin = value.sd_below, ymax = value.sd_above), 
                alpha = .4, color = 'black') + 
    geom_line(aes(d, value.mean), alpha = .9, color = 'black')  + 
    geom_vline(aes(xintercept = xshift), 
               color = 'grey50', linetype = 'dashed',
               data = shift_df) + 
    facet_grid(indicator ~ pretty_delta, scales = 'free_y', 
               switch = "y") + 
    ggtitle('Forestgap model') + 
    xlab( expression(m[0]) ) + 
    theme_ifcam
    
ggsave(paste0(output_figure_path, FILE_PREFIX, "forestgap.pdf"),
       width = GRAPH_WIDTH, height = GRAPH_HEIGHT,
       plot = forestgap_csplot)

rm(branch)
gc()







# Musselbed model
# -------------------------

if ( !exists("branch") ) { 
  load_url(files[['musselbed_cs']], verbose = TRUE)
}

upper_result <- add_indics(branch[['DatBif']], 
                           branch[['snaps']], 
                           branch[['time_series']],
                           subsize = CG_SUBSIZE, 
                           null_replicates = NULL_REPS, 
                           .progress = 'time',
                           .parallel = FALSE)
# save.image("musselbed.dump.rda")

upper_result_trend <- summarise_trend(upper_result, ~ ID + indicator + delta + d)

shift_df <- create_shift_data_frame(upper_result_trend, "d", "delta")

upper_result_trend[ ,'pretty_d'] <- with(upper_result_trend, paste0("d = ", d))
shift_df[ ,'pretty_d'] <- with(shift_df, paste0("d = ", d))

save(upper_result, upper_result_trend, shift_df, 
     file = paste0(output_figure_path, "musselbed_all_indic_", 
                   CG_SUBSIZE, "x", CG_SUBSIZE, ".rda"))
     

# save.image('grazing_allindic.rda')

# Plot chosen cross-sections
musselbed_csplot <- 
  ggplot(upper_result_trend) + 
    geom_ribbon(aes(delta, ymin = null.sd_below, ymax = null.sd_above), 
                alpha = .2, color = 'red') + 
    geom_line(aes(delta, null.mean), color = 'red') + 
    geom_ribbon(aes(delta, ymin = value.sd_below, ymax = value.sd_above), 
                alpha = .4, color = 'black') + 
    geom_line(aes(delta, value.mean), alpha = .9, color = 'black')  + 
    geom_vline(aes(xintercept = xshift), 
               color = 'grey50', linetype = 'dashed',
               data = shift_df) + 
    facet_grid(indicator ~ pretty_d, scales = 'free_y', 
               switch = "y") + 
    ggtitle('Musselbed model') + 
    xlab( expression(m[0]) ) + 
    theme_ifcam

ggsave(paste0(output_figure_path, FILE_PREFIX, "musselbed.pdf"),
       width = GRAPH_WIDTH, height = GRAPH_HEIGHT,
       plot = musselbed_csplot)

rm(branch)
gc()

