
# 
# 

# Setwd
setwd('/home/alex/work/2014-2015/SpatialStress/ifcam_output/figures/generic_indicators/grazing_zoom')
# setwd('/home/alex/system/tmp')

library(caspr)
library(plyr)
library(spatialwarnings)
library(ggplot2)


# 
# # Use this to use the proxy at IISC
# Sys.setenv(https_proxy = "proxy.iisc.ernet.in:3128")
# Sys.setenv(http_proxy  = "proxy.iisc.ernet.in:3128")

# TODO: 
#   - add time series (done)
#   - set t_min to 500 (done)
#   - set t_max to 3000 (done)
#   - forestgap: d btw 0 and 0.25, delta between 0a nd 0.25 (done)
#     * NOTE: no, I put delta to 0.35 as it allows showing the whole 
#         trends. 
#   - musselbed: removed the value \delta = 0 from plots (nothing to do here)
#   - musselbed: 0 to 0.3 (done)
#   - musselbed: catastrophic shifts along the b parameter ? (nope)
#   - grazing: use g and m (done)
# 

# Load package/install if necessary
library(devtools)
install_github('fdschneider/caspr') # might need a restart of R
library(caspr)
library(ggplot2)
library(doParallel)

# Global variables
# --------------------------------------------------------
SIZE_CS <- 400  # Size of lattice for cross-sections
RES_CS  <- 20  # Number of points on each dimension
NSNAPS  <- 30   # Number of snapshots to save at the end of simulation
LENGTH.STAT <- 200 # Number of time steps to consider when computing mean covers
TMIN <- 30000  # Minimum simulation time before taking snaps 
TMAX <- 40000 # Maximum simulation time before taking snaps (if no eq is found)

# Computation-related variables
PARALLEL <- TRUE
NCORES   <- 5


# Model-specific variables
# --------------------------------------------------------
# Initial covers 
GRAZING_INIT_UPPER <- c(.99, .01, 0) # veg/barren/disturbed (+/0/-)
GRAZING_INIT_LOWER <- c(.05, .95, 0)
FORESTGAP_INIT_UPPER <- c(.95, .05) # veg/gap 
FORESTGAP_INIT_LOWER <- c(.05, .95)
MUSSELBED_INIT_UPPER <- c(.99, .01, 0) # veg/gap/disturbed
MUSSELBED_INIT_LOWER <- c(.05, .95, 0)


# Parameters
GRAZING_DEFAULT_PARMS <- list(del = 0.9, b = 0.5, c_ = 0.2, r = 0.01, 
                              f = 0.9, d = 0.1, p   = 1, 
                              # changed to sequences later on
                              m0  = .2, # homogeneous
                              g   = .3) # spatial

GRAZING_CS_PARMS   <- list(m0 = seq(0.08, .1, length.out = RES_CS),
                           g  = c(0.18))



# Helper function that runs simulation, taking into account global params
do_simus <- function(model, parms, init, size) { 
  
  # We reset the cluster with each call to do_simus, otherwise 
  #   memory consumptions explodes
  # Register parallel backend
  if (PARALLEL) { 
    stopImplicitCluster()
    registerDoParallel(cores = NCORES)
  } else { 
    registerDoSEQ()
  }
  
  ca_arraySS(model, init = init, width = size, height = size, 
             parms = parms, nsnaps = NSNAPS, length.stat = LENGTH.STAT,
             t_min = TMIN, t_max = TMAX)
}




# Grazing model -------------
parms <- GRAZING_DEFAULT_PARMS
parms[names(GRAZING_CS_PARMS)] <- GRAZING_CS_PARMS # mind the single bracket

result_grazing_upper <- do_simus(grazing, parms, GRAZING_INIT_UPPER, SIZE_CS)

save(result, 
     file = "./result_grazing_zoom.rda", compress = 'bzip2')

gc()


# Compute spatial indicators
# Extract matrices
matrices <- llply(result[[2]], llply, caspr:::as.matrix.landscape)
matrices <- llply(matrices, llply, function(x) x == "+") # keep only + state

# Compute indicators
genics <- llply(matrices, generic_spews, 
                detrend = FALSE, subsize = 4, moranI_coarse_grain = FALSE, 
                .progress = "time")

genitest <- llply(genics, indictest, 
                  null_replicates = 9, 
                  .progress = 'time')

plotdat <- llply(genitest, as.data.frame)
plotdat <- Map(function(n,e) data.frame(ID = n, e), 
               seq.int(length(plotdat)), plotdat)
plotdat <- do.call(rbind, plotdat)
plotdat <- join(result[["DatBif"]], plotdat, by = "ID")

# Summarize trend
plotdat_trend <- ddply(plotdat, ~ ID + indicator + g + m0, 
                       summarise, 
                       value.mean = mean(value), 
                       value.sd_above = mean(value) + sd(value), 
                       value.sd_below = mean(value) - sd(value),
                       null.mean = mean(null_mean), 
                       null.sd_above = mean(null_mean) + sd(null_mean), 
                       null.sd_below = mean(null_mean) - sd(null_mean),
                       .progress = 'time')



ggplot(plotdat_trend) + 
  geom_ribbon(aes(m0, ymin = null.sd_below, ymax = null.sd_above), alpha = .2, 
              color = 'red') + 
  geom_line(aes(m0, null.mean), color = 'red') + 
  geom_ribbon(aes(m0, ymin = value.sd_below, ymax = value.sd_above), 
              alpha = .4, color = 'black') + 
  geom_line(aes(m0, value.mean), alpha = .6, color = 'black')  + 
  facet_grid( indicator ~ g, scales = 'free_y', labeller = label_both) + 
  ggtitle('Musselbed model')+
  annotate(geom = "vline", xintercept = .089, x = .089)

ggsave('./grazing_genic_zoom_no_null.pdf', width = 6, height = 10)


  


