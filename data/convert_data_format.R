# 
# This function converts raw data from the output of ca_snapsSS et al. to one 
#   that the spatialwarnings can understand. 
# 

# Base paths (with / at the end /!\)
result_folder <- '/home/alex/work/2014-2015/SpatialStress/ifcam_output/data/' 
setwd(result_folder)

# We use the binary_matrix as defined in caspr as it provides a landscape method
library(caspr)

CONVERT_PD <- TRUE
CONVERT_CS <- TRUE

# Function that actually does the conversion
convert_dat <- function(upper_branch, lower_branch, newfile, state) { 

  # Convert all matrices to binary_matrices (for loop = less memory used)
  for (i in seq.int(upper_branch[[2]])) { 
    cat('.')
    upper_branch[[2]][[i]] <- lapply(upper_branch[[2]][[i]], 
                                     as.binary_matrix, is = state)
    lower_branch[[2]][[i]] <- lapply(lower_branch[[2]][[i]], 
                                     as.binary_matrix, is = state)
  }
  cat('\n')
  save(upper_branch, lower_branch, file = newfile, compress = 'bzip2')
  
  rm(lower_branch, upper_branch)
}

# Function that actually does the conversion, just for a single branch
convert_dat_onebranch <- function(branch, newfile, state) { 

  # Convert all matrices to binary_matrices (for loop = less memory used)
  for (i in seq.int(branch[[2]])) { 
    cat('.')
    branch[[2]][[i]] <- lapply(branch[[2]][[i]], as.binary_matrix, is = state)
  }
  cat('\n')
  save(branch, file = newfile, compress = 'bzip2')
  
  rm(branch)
}






# Convert phase-diagrams data
# ------------------------------------

if ( CONVERT_PD ) { 

  # Convert forestgap 
  load(paste0(result_folder,'result_forestgap.rda'), verbose = TRUE)
  convert_dat(result_forestgap_upper, result_forestgap_lower, 
              newfile = paste0(result_folder, "result_forestgap_processed.rda" ),
              state = '+')
  rm(result_forestgap_upper, result_forestgap_lower)
  gc()
  
  # Convert musselbed
  load(paste0(result_folder,'result_musselbed.rda'), verbose = TRUE)
  convert_dat(result_musselbed_upper, result_musselbed_lower, 
              newfile = paste0(result_folder, "result_musselbed_processed.rda" ),
              state = '+')
  rm(result_musselbed_upper, result_musselbed_lower)
  gc()
  
  # Convert grazing
  load(paste0(result_folder,'result_grazing.rda'), verbose = TRUE)
  convert_dat(result_grazing_upper, result_grazing_lower, 
              newfile = paste0(result_folder, "result_grazing_processed.rda" ),
              state = '+')
  rm(result_grazing_upper, result_grazing_lower)
  gc()
  
}


# Convert cross-section data
# ------------------------------------

if ( CONVERT_CS ) { 

  # Convert musselbed
  load(paste0(result_folder,'result_musselbed_cs.rda'), verbose = TRUE)
  convert_dat_onebranch(result_musselbed_upper, 
                        newfile = paste0(result_folder, "result_musselbed_cs_processed.rda" ), 
                        state = '+')
  rm(result_musselbed_upper)
  gc()
  
  # Convert forestgap
  load(paste0(result_folder,'result_forestgap_cs.rda'), verbose = TRUE)
  convert_dat_onebranch(result_forestgap_upper, 
                        newfile = paste0(result_folder, "result_forestgap_cs_processed.rda" ), 
                        state = '+')
  rm(result_forestgap_upper)
  gc()
  
  # Convert grazing
  load(paste0(result_folder,'result_grazing_cs.rda'), verbose = TRUE)
  convert_dat_onebranch(result_grazing_upper, 
                        newfile = paste0(result_folder, "result_grazing_cs_processed.rda" ), 
                        state = '+')
  rm(result_grazing_upper) 
  gc()

}
