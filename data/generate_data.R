# 
# Generate the data needed for the analyses 

# Setwd
setwd('~/caspr/draft/ifcam_results')
# setwd('/home/alex/system/tmp')

# 
# # Use this to use the proxy at IISC
# Sys.setenv(https_proxy = "proxy.iisc.ernet.in:3128")
# Sys.setenv(http_proxy  = "proxy.iisc.ernet.in:3128")


# Load package/install if necessary
library(devtools)
# install_github('fdschneider/caspr') # might need a restart of R
library(caspr)
library(ggplot2)
library(doParallel)

# Global variables
SIZE_PD    <- 100 # Size of lattice for 2D plane diagram
SIZE_CS <- 400    # Size of lattice for cross-sections
RES     <- 51     # Number of points on each dimension
RES_CS  <- 201    # Number of points on each dimension
NSNAPS <- 10      # Number of snapshots to save at the end of simulation
REDO_COMPUTATIONS_PD <- FALSE # Redo computations for phase diagrams
REDO_COMPUTATIONS_CS <- TRUE # Redo computations for cross-sections
LENGTH.STAT <- 200 # Number of time steps to consider when computing mean covers
TMIN <- 500  # Minimum simulation time
TMAX <- 3000 # Maximum simulation time (in case no eq is found)

# Computation-related variables
PARALLEL <- TRUE
NCORES   <- 23

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

# Register parallel backend
if (PARALLEL) { 
  registerDoParallel(cores = NCORES)
} else { 
  registerDoSEQ()
}

# registerDoSEQ()

# Produce data needed for the 2D state diagram
# --------------------------------------------------------
if (REDO_COMPUTATIONS_PD) {
  
  
  # Grazing model -------------
  parms <- list(del = 0.9,
                b   = 0.5, # suggested default value in the model
                c_  = 0.2,
                m0  = seq(0, .4, length.out = RES),
                g   = seq(0, .6, length.out = RES),
                r   = 0.01,
                f   = 0.9,
                d   = 0.1,
                p   = 1)
  
  result_grazing_upper <- ca_arraySS(grazing, 
                                     init = c(.99, .01, 0), 
                                     width = SIZE_PD, height = SIZE_PD, 
                                     parms = parms, nsnaps = NSNAPS,
                                     length.stat = LENGTH.STAT,
                                     t_min = TMIN, t_max = TMAX)
  
  result_grazing_lower <- ca_arraySS(grazing, 
                                     init = c(.05, .95, 0), 
                                     width = SIZE_PD, height = SIZE_PD, 
                                     parms = parms, nsnaps = NSNAPS,
                                     length.stat = LENGTH.STAT,
                                     t_min = 500, t_max = 3000,
                                     t_min = TMIN, t_max = TMAX)
  
  save(result_grazing_lower, result_grazing_upper, 
       file = "./result_grazing.rda", compress = 'bzip2')
  rm(result_grazing_upper)
  rm(result_grazing_lower)
  
  
  # Forestgap model -------------
  parms <- list(alpha = 0.2,
                delta = seq(0, .35, length.out = RES),
                d     = seq(0, .25, length.out = RES))
  
  result_forestgap_upper <- ca_arraySS(forestgap, 
                                       init = c(.95, .05), 
                                       width = SIZE_PD, height = SIZE_PD, 
                                       parms = parms, nsnaps = NSNAPS,
                                       length.stat = LENGTH.STAT,
                                       t_min = TMIN, t_max = TMAX)
  
  result_forestgap_lower <- ca_arraySS(forestgap, 
                                       init = c(.05, .95), 
                                       width = SIZE_PD, height = SIZE_PD, 
                                       parms = parms, nsnaps = NSNAPS,
                                       length.stat = LENGTH.STAT,
                                       t_min = TMIN, t_max = TMAX)
  
  save(result_forestgap_lower, result_forestgap_upper, 
       file = "./result_forestgap.rda", compress = 'bzip2')
  
  rm(result_forestgap_lower)
  rm(result_forestgap_upper)
  
  
  # Musselbed model -------------
  parms <- list(r     = 0.4,
                d     = seq(0, 1,  length.out = RES),
                delta = seq(0, .3, length.out = RES))
  
  result_musselbed_upper <- ca_arraySS(musselbed, 
                                       init = c(.99, .01, 0), 
                                       width = SIZE_PD, height = SIZE_PD, 
                                       parms = parms, nsnaps = NSNAPS,
                                       length.stat = LENGTH.STAT,
                                       t_min = TMIN, t_max = TMAX)
  
  result_musselbed_lower <- ca_arraySS(musselbed, 
                                       init = c(.05, .95, 0), 
                                       width = SIZE_PD, height = SIZE_PD, 
                                       parms = parms, nsnaps = NSNAPS,
                                       length.stat = LENGTH.STAT,
                                       t_min = TMIN, t_max = TMAX)
  
  save(result_musselbed_lower, result_musselbed_upper, 
       file = "./result_musselbed.rda", compress = 'bzip2')
  
  rm(result_musselbed_lower)
  rm(result_musselbed_upper)
    
} 

# Produce data needed for cross-sections
# --------------------------------------------------------
if (REDO_COMPUTATIONS_CS) { 
  
  
  
  
}


