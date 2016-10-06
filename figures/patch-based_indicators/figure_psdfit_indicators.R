# 
# 
# This file contains code that plots spectral indicators along the stress 
#   gradient. Spectral indicators include the r-spectrum and SDR that 
#   derives from it.
# 

library(ggplot2)
library(plyr) 
library(tidyr) 
library(grid) 
library(gridExtra) # for grid.arrange
library(doParallel) 
library(devtools) 

# Install latest spatialwarnings
# install_github('fdschneider/spatial_warnings') # mind the _
install_github('fdschneider/spatial_warnings') # mind the _
library(spatialwarnings) 

# Define output path 
# /!\ NB: folders must have a trailing /
working_directory  <- '/home/alex/work/2014-2015/SpatialStress/ifcam_output/figures/'
setwd(working_directory)


# Graphic output 
GRAPH_WIDTH <- 9
GRAPH_HEIGHT <- 3  # The 2x value ensures that 2D planes are squareish
FILE_PREFIX  <- paste0('psdtype_indicator_')
output_figure_path <- paste0(working_directory, 'patch-based_indicators/')

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
  registerDoParallel(cores = 16)
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
add_psdfit <- function(df, # data.frame of all simulations (DatBif)
                       matrices, # corresponding matrices (list of lists)
                       ...) { 
  
  # Get psd types (this is brittle)
  result <- llply(matrices, indicator_psdtype, merge = TRUE, 
                  .progress = 'time', 
                  .parallel = PARALLEL,  
                  ...)
  
  # Map with ids
  result <- Map(function(n, e) data.frame(ID = n, e), 
                seq.int(length(result)), result)
  result <- do.call(rbind, result)
  
  # Join with datbif 
  result <- join(df, result, by = "ID")
  
  
  return(result)
}



# Grazing model
# -------------------------
if ( ! exists('upper_branch') ) { 
  load_url(files[['grazing']], verbose = TRUE)
  rm('lower_branch'); gc()
}


result <- add_psdfit(upper_branch[['DatBif']],
                     upper_branch[['snaps']])

save.image("psds.gr.dump.rda")
# load("psds.gr.dump.rda")

ggplot() + 
  geom_tile(aes(x = m0, y = g, fill = type), 
            data = subset(result, best)) + 
  geom_tile(aes(x = m0, y = g), fill = "gray50", 
            data = subset(result, best & percolation == 1))

# Remove lognormal
result.noLnorm <- subset(result, type != 'lnorm')

result.all <- ddply(result.noLnorm, ~ m0 + g, function(df) { 
#                     if (!is.na(df$method)) browser()
                    return(with(df, 
                      data.frame(aic = type[AIC == min(AIC)],
                                 aicc = type[AICc == min(AICc)],
                                 bic = type[BIC == min(BIC)])))
                    })

result.all <- gather(result.all, crit, type, aic, aicc, bic)

ggplot() + 
  geom_tile(aes(x = m0, y = g, fill = type), 
            data = subset(result.all)) + 
  geom_tile(aes(x = m0, y = g), fill = "gray50", 
            data = subset(result, best & percolation == 1)) + 
  facet_wrap( ~ crit) + 
  scale_fill_manual(values = c('#E8643A','#34AA35','#E4CA1B','#26A9E2'))
  
ggsave(paste0(output_figure_path, "grazing_psdtype_draft_nolnorm.pdf"), 
       width = 10, height = 4)




# Forestgap model
# -------------------------
if ( ! exists('upper_branch') ) { 
  load_url(files[['forestgap']], verbose = TRUE)
  rm('lower_branch'); gc()
}


result <- add_psdfit(upper_branch[['DatBif']],
                     upper_branch[['snaps']])

save.image("psds.fg.dump.rda")
# load("psds.fg.dump.rda")

# Remove lognormal
result.noLnorm <- subset(result, type != 'lnorm')

result.all <- ddply(result.noLnorm, ~ d + delta, function(df) { 
                    return(with(df, 
                      data.frame(aic = type[AIC == min(AIC)],
                                 aicc = type[AICc == min(AICc)],
                                 bic = type[BIC == min(BIC)])))
                    })

result.all <- gather(result.all, crit, type, aic, aicc, bic)

ggplot() + 
  geom_tile(aes(x = d, y = delta, fill = type), 
            data = subset(result.all)) + 
  geom_tile(aes(x = d, y = delta), fill = "gray50", 
            data = subset(result, best & percolation == 1)) + 
  facet_wrap( ~ crit) + 
  scale_fill_manual(values = c('#E8643A','#34AA35','#E4CA1B','#26A9E2'))
  
ggsave(paste0(output_figure_path, "forestgap_psdtype_draft_nolnorm.pdf"), 
       width = 10, height = 4)




# Musselbed model
# -------------------------
if ( ! exists('upper_branch') ) { 
  load_url(files[['musselbed']], verbose = TRUE)
  rm('lower_branch'); gc()
}


result <- add_psdfit(upper_branch[['DatBif']],
                     upper_branch[['snaps']])

save.image("psds.mb.dump.rda")
# load("psds.mb.dump.rda")

# Remove lognormal
result.noLnorm <- subset(result, type != 'lnorm')

result.all <- ddply(result.noLnorm, ~ delta + d, function(df) { 
                    return(with(df, 
                      data.frame(aic = type[AIC == min(AIC)],
                                 aicc = type[AICc == min(AICc)],
                                 bic = type[BIC == min(BIC)])))
                    })

result.all <- gather(result.all, crit, type, aic, aicc, bic)

ggplot() + 
  geom_tile(aes(x = delta, y = d, fill = type), 
            data = subset(result.all)) + 
  geom_tile(aes(x = delta, y = d), fill = "gray50", 
            data = subset(result, best & percolation == 1)) + 
  facet_wrap( ~ crit) + 
  scale_fill_manual(values = c('#E8643A','#34AA35','#E4CA1B','#26A9E2'))
  
ggsave(paste0(output_figure_path, "musselbed_psdtype_draft_nolnorm.pdf"), 
       width = 10, height = 4)

