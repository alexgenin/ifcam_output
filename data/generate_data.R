# 
# Generate the data needed for the analyses 

# Setwd
setwd('/home/alex/work/2014-2015/SpatialStress/ifcam_output/data/')
# setwd('/home/alex/system/tmp')

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
# install_github('fdschneider/caspr') # might need a restart of R
library(caspr)
library(ggplot2)
library(doParallel)

# Global variables
# --------------------------------------------------------
SIZE_PD <- 100  # Size of lattice for 2D plane diagram
SIZE_CS <- 400  # Size of lattice for cross-sections
RES_PD  <- 51   # Number of points on each dimension
RES_CS  <- 201  # Number of points on each dimension
NSNAPS  <- 10   # Number of snapshots to save at the end of simulation
REDO_COMPUTATIONS_PD <- FALSE # Redo computations for phase diagrams
REDO_COMPUTATIONS_CS <- TRUE # Redo computations for cross-sections
LENGTH.STAT <- 200 # Number of time steps to consider when computing mean covers
TMIN <- 500  # Minimum simulation time
TMAX <- 3000 # Maximum simulation time (in case no eq is found)

# Computation-related variables
PARALLEL <- TRUE
NCORES   <- 23


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

FORESTGAP_DEFAULT_PARMS <- list(alpha = 0.2, 
                                # changed to sequences later on 
                                delta = 0.17, d = 0.125)

MUSSELBED_DEFAULT_PARMS <- list(r = 0.4, 
                                # changed to sequences later on
                                d = 0.5, delta = 0.15)


# 2D Phase diagram parameters
GRAZING_PD_PARMS   <- list(m0 = seq(0, .4, length.out = RES_PD), 
                           g =  seq(0, .6, length.out = RES_PD))

MUSSELBED_PD_PARMS <- list(d =     seq(0, 1,  length.out = RES_PD), 
                           delta = seq(0, .3, length.out = RES_PD))

FORESTGAP_PD_PARMS <- list(delta = seq(0, .35, length.out = RES_PD), 
                           d =     seq(0, .25, length.out = RES_PD))


# Cross-sections parameters
GRAZING_CS_PARMS   <- list(m0 = seq(0, .4, length.out = RES_CS),
                           g  = c(0.18, 0.072, 0))

MUSSELBED_CS_PARMS <- list(d = c(0.60, 0.24, 0), 
                           delta = seq(0, .3, length.out = RES_PD))

FORESTGAP_CS_PARMS <- list(delta = c(0.189, 0.07, 0), 
                           d = seq(0, .25, length.out = RES_PD))


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





# registerDoSEQ()

# Produce data needed for the 2D state diagram
# --------------------------------------------------------
if (REDO_COMPUTATIONS_PD) {
  
  
  
  # Grazing model -------------
  parms <- GRAZING_DEFAULT_PARMS
  parms[names(GRAZING_PD_PARMS)] <- GRAZING_PD_PARMS # mind the single bracket
  
  result_grazing_upper <- do_simus(grazing, parms, GRAZING_INIT_UPPER, SIZE_PD)
  result_grazing_lower <- do_simus(grazing, parms, GRAZING_INIT_LOWER, SIZE_PD)
    
  save(result_grazing_lower, result_grazing_upper, 
       file = "./result_grazing.rda", compress = 'bzip2')
  rm(result_grazing_upper)
  rm(result_grazing_lower)
  
  
  
  # Forestgap model -------------
  parms <- FORESTGAP_DEFAULT_PARMS
  parms[names(FORESTGAP_PD_PARMS)] <- FORESTGAP_PD_PARMS

  result_forestgap_upper <- do_simus(grazing, parms, FORESTGAP_INIT_UPPER, SIZE_PD)
#   result_forestgap_lower <- do_simus(grazing, parms, FORESTGAP_INIT_LOWER, SIZE_PD)
  
  save(result_forestgap_lower, result_forestgap_upper, 
       file = "./result_forestgap.rda", compress = 'bzip2')
  
  rm(result_forestgap_lower)
  rm(result_forestgap_upper)
  
  
  
  # Musselbed model -------------
  
  parms <- MUSSELBED_DEFAULT_PARMS
  parms[names(MUSSELBED_PD_PARMS)] <- MUSSELBED_PD_PARMS
  
  result_musselbed_upper <- do_simus(grazing, parms, MUSSELBED_INIT_UPPER, SIZE_PD)
#   result_musselbed_lower <- do_simus(grazing, parms, MUSSELBED_INIT_LOWER, SIZE_PD)
  
  save(result_musselbed_lower, result_musselbed_upper, 
       file = "./result_musselbed.rda", compress = 'bzip2')
  
  rm(result_musselbed_lower)
  rm(result_musselbed_upper)
    
} 




# Produce data needed for cross-sections
# --------------------------------------------------------

if (REDO_COMPUTATIONS_CS) { 
  
  
  # Grazing model -------------
  parms <- GRAZING_DEFAULT_PARMS
  parms[names(GRAZING_CS_PARMS)] <- GRAZING_CS_PARMS # mind the single bracket
  
#   result_grazing_upper <- do_simus(grazing, parms, GRAZING_INIT_UPPER, SIZE_CS)
  result_grazing_lower <- do_simus(grazing, parms, GRAZING_INIT_LOWER, SIZE_CS)
  
  save(result_grazing_lower, result_grazing_upper, 
       file = "./result_grazing_cs.rda", compress = 'bzip2')
  
  rm(result_grazing_upper)
  rm(result_grazing_lower)
  gc()
  
  # Forestgap model -------------
  parms <- FORESTGAP_DEFAULT_PARMS
  parms[names(FORESTGAP_CS_PARMS)] <- FORESTGAP_CS_PARMS 
  
  result_forestgap_upper <- do_simus(forestgap, parms, FORESTGAP_INIT_UPPER, SIZE_CS)
#   result_forestgap_lower <- do_simus(forestgap, parms, FORESTGAP_INIT_LOWER, SIZE_CS)
  
  save(result_forestgap_upper, 
      file = "./result_forestgap_cs.rda", compress = 'bzip2')
  
#   rm(result_forestgap_lower)
  rm(result_forestgap_upper)
  gc()
  
  # Musselbed model -------------
  parms <- MUSSELBED_DEFAULT_PARMS
  parms[names(MUSSELBED_CS_PARMS)] <- MUSSELBED_CS_PARMS 
  
  result_musselbed_upper <- do_simus(musselbed, parms, MUSSELBED_INIT_UPPER, SIZE_CS)
  result_musselbed_lower <- do_simus(musselbed, parms, MUSSELBED_INIT_LOWER, SIZE_CS)
  
  save(result_musselbed_lower, result_musselbed_upper, 
      file = "./result_musselbed_cs.rda", compress = 'bzip2')
  
  rm(result_musselbed_lower)
  rm(result_musselbed_upper)
  gc()
  
}


