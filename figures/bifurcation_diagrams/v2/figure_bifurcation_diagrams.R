# 
# 
# This file contains code that generate bifurcation diagrams for all models. 
# 
# It generates 2D planes with contours for the upper branch (see Flo's 
#   manuscript, Fig 2.) and cross-sections chose for certain thresholds.
# 

# Some prerequisite before anything else
# Set proxy if needed

# Load libraries
library(ggplot2)
library(grid) 
library(gridExtra) # for grid.arrange
library(plyr) 
library(rootSolve)
library(tidyr)

# Define output path 
# /!\ NB: folders must have a trailing /
working_directory  <- '/home/alex/work/2014-2015/SpatialStress/ifcam_output/figures/'
output_figure_path <- "/home/alex/work/2014-2015/SpatialStress/ifcam_output/figures/bifurcation_diagrams/v2/" 

# Define data paths 
# /!\ NB: folders must have a trailing /
data_folder <- '/home/alex/work/2014-2015/SpatialStress/ifcam_output/data/'
# Use this folder if you do not have the results on your computer: it will fetch
#   the latest ones.
# data_folder <- "http://alex.lecairn.org/ifcam/" 
files <- list(musselbed_cs = paste0(data_folder, "result_musselbed_cs_processed.rda"),
              grazing   = paste0(data_folder, "result_grazing_processed.rda"),
              forestgap = paste0(data_folder, "result_forestgap_processed.rda"))




# Define some helper functions

# Graph theme
theme_ifcam <- 
  theme_bw() + 
  theme(strip.background = element_blank(), 
        panel.border = element_rect(fill = NA, linetype = 'dotted', color = 'grey20'))
  
no_legend <- theme(legend.position = "none")


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

# Rank values in a different way than rank
rank2 <- .%>% rank %>% as.factor %>% as.integer %>% as.character

# Merge two summary data.frames
merge_branches <- function(upper, lower) { 
  data.frame(branch = c(rep('Upper branch', nrow(upper)), 
                        rep('Lower branch', nrow(lower))),
             rbind(upper, lower))
}

# Get the closest value in a vector
closest_to <- function(val, X, quiet = FALSE) { 
  index <- which( abs(val - X) == min(abs(val - X)) )
  new_value <- X[min(index)]
  if (!quiet) { 
    cat('Picked value ', new_value, ' (error: ', new_value - val,')\n', sep ='')
  }
  return(new_value)
}


# Merge the plots in a three panes layout
merge_plots <- function(phasediag, cs_high, cs_mid, cs_low) { 
  grid.arrange(phasediag, 
               cs_high, cs_mid, cs_low, 
               ncol = 1, 
               heights = c(.5, .5/3, .5/3, .5/3))
}

# Values common for all models 
setwd(working_directory)

# Graphic output 
GRAPH_WIDTH <- 5
GRAPH_HEIGHT <- GRAPH_WIDTH * 2 # The 2x value ensures that 2D planes are squareish
FILE_PREFIX <- 'bifurc_diagram_'






# Grazing model 
# ------------------


# Load data regarding spatial model
if ( ! exists("upper_branch") ) { load_url(files[["grazing"]]) }

datbif <- merge_branches(upper_branch[["DatBif"]], 
                         lower_branch[["DatBif"]])
rm(upper_branch, lower_branch)

# Subset datbig with the values we use for spatial simulations
g_values <- sapply(c(0, 0.072, 0.18), closest_to, datbif[ ,'g'])
datbif_spatial <- subset(datbif, g %in% g_values)

datbif_spatial[ ,'rhop'] <- datbif_spatial[ ,'mean_cover_.'] # give it a sensible nameà)
datbif_spatial[ ,'eq_type'] <- "stable" # we only have stable eqs in spatial sims
datbif_spatial[ ,'hetero'] <- rank2(datbif_spatial[ ,'g'])
datbif_spatial[ ,'hetero_rank'] <- rank2(datbif_spatial[ ,'g'])



# Find roots using mean-field models equations
# Common grazing parameters used
r <- 0.01; f <- .9; b <- .5; c <- .2; d <- .5; # d <- .1 in actual runs
grazing_zpoly <- function(g, m0) { # coefficient-generating function
  c( (b*r - m0*d - m0*r - g*d - g*r) / c * f,  # degree 0
      b/c - m0/c - g/c - b*r/(c*f) - r/f + (g*d+g*r)/c*f, # degree 1
      r/f - b/c + g/c - 1, # degree 2
      1) # degree 3
}
  
# Given a polynomial coefficients, find its roots and format them 
#   in a nice way for plotting
getroots_grazing <- function(zpoly, pars) { 
  ddply(pars, names(pars), function(df) { 
      z <- zpoly(df[ ,1], df[ ,2])
      roots <- Re(polyroot(z))
      
      # We discard two first roots if two of them are repeated (with a tol)
      if ( length(roots) >= 3 && abs(roots[1] - roots[2]) < 1e-10) { 
        roots <- roots[3]
      }
      
      # Identify the type of stability
      stability <- ifelse(length(roots) > 2, "bistable", "monostable")
      
      # Identiy eq type
      if ( length(roots) > 2 ) { 
        eq_type <- c("unstable", "stable", "unstable")
        branch <- c('lower', "unst", 'upper')
      } else { 
        eq_type <- c('unstable', "unstable")
        branch  <- c("unst", 'upper')
      }
      
      data.frame(branch = branch, 
                 stability = stability, eq_type = eq_type,
                 rhop = roots)
  })
}

# Use g values for 
graph_values <- expand.grid(g = c(0, 0.1, 0.18),
                            m0 = seq(0, .4, length.out = 500))
datbif_meanf <- getroots_grazing(grazing_zpoly, graph_values)

datbif_meanf[ ,'hetero_rank'] <- rank2(datbif_meanf[ ,'g'])


# Merge both datasets
datbif_grazing <- rbind.fill(data.frame(type = "Mean-field approx.", datbif_meanf), 
                             data.frame(type = "Spatial simulation", 
                                        eq_type = "stable", # we only have stable eqs in 
                                        datbif_spatial))

datbif_grazing[ ,'hetero'] <- datbif_grazing[ ,'g']
datbif_grazing[ ,'homo'] <- datbif_grazing[ ,'m0']
datbif_grazing[ ,'pretty_hetero'] <- with(datbif_grazing, paste0("g = ", g))
datbif_grazing[ ,'model'] <- "Grazing model"







# Forestgap model 
# ------------------


# Load data regarding spatial model
if ( ! exists("upper_branch") ) { load_url(files[["forestgap"]]) }

datbif <- merge_branches(upper_branch[["DatBif"]], 
                         lower_branch[["DatBif"]])
rm(upper_branch, lower_branch)

delta_values <- sapply(c(0, 0.07, 0.189), closest_to, datbif[ ,'delta'])
datbif_spatial <- subset(datbif, delta %in% delta_values)
datbif_spatial[ ,'rhop'] <- datbif_spatial[ ,'mean_cover_.'] # give it a sensible nameà)
datbif_spatial[ ,'eq_type'] <- "stable" # we only have stable eqs in spatial sims
datbif_spatial[ ,'hetero_rank'] <- rank2(datbif_spatial[ ,'delta'])



# Common forestgap parameters used
alpha <- 1; # b in Sabiha's document ?
graph_values <- expand.grid(delta = c(0.5, 1, 1.5),
                            d = seq(0, 1, length.out = 500))

# Note: we have an explicit solution for \rho+ and \rho-, so we juste use that
# Note2: there is a typo in the determinant in Sabiha's doc
datbif_meanf <- ddply(graph_values, names(graph_values), function(df) { 
  delta <- df[ ,'delta']; d <- df[ ,'d']
  det <- abs(sqrt( (alpha + delta + d)^2 - 4*alpha*delta )) # typo in doc here
  rbind( data.frame(eq_type = "unstable", # see document
                    branch = NA, 
                    rhop = ((alpha + delta + d) + det) / (2 * alpha)), 
         data.frame(eq_type = "stable", # see document
                    branch = NA, 
                    rhop = ((alpha + delta + d) - det) / (2 * alpha)) )
})
         
datbif_meanf[ ,'hetero_rank'] <- rank2(datbif_meanf[ ,'delta'])

# Merge both datasets
datbif_forestgap <- rbind.fill(data.frame(type = "Mean-field approx.", datbif_meanf), 
                     data.frame(type = "Spatial simulation", 
                                eq_type = "stable", # we only have stable eqs in 
                                datbif_spatial))

datbif_forestgap[ ,'homo'] <- datbif_forestgap[ ,'d']
datbif_forestgap[ ,'hetero'] <- datbif_forestgap[ ,'delta']
datbif_forestgap[ ,'pretty_hetero'] <- with(datbif_forestgap, paste0("delta = ", delta))
datbif_forestgap[ ,'model'] <- "Forest Gap Model"






# Musselbed model 
# ------------------


# Load data regarding spatial model
if ( ! exists("upper_branch") ) { load_url(files[["musselbed_cs"]]) }

dat <- branch[["DatBif"]]

# Note that we need .6 here but we only have up to .3 -> rerun simulations
d_values <- sapply(c(0, .24, .6), closest_to, dat[ ,'d'])
datbif_spatial <- subset(dat, d %in% d_values)
datbif_spatial[ ,'rhop'] <- datbif_spatial[ ,'mean_cover_.'] # give it a sensible nameà)
datbif_spatial[ ,'eq_type'] <- "stable" # we only have stable eqs in spatial sims
datbif_spatial[ ,'hetero_rank'] <- rank2(datbif_spatial[ ,'d'])


# Common musselbed parameters used
alpha2 <- 1 # b in Sabiha's document ?alpha
graph_values <- expand.grid(d = c(0, .075, .15, .4), 
                            delta = seq(0, 1, length.out = 100))

# We reuse the values used in the mean field doc
datafiles <- dir('./bifurcation_diagrams/v2/mussel_dat/', full = TRUE)

datbif_meanf <- alply(datafiles, 1, function(name) { 
       dat <- read.table(name)
       names(dat) <- c('delta', 'root1', 'root2')
       gather(dat, root, rhop, root1, root2)
     })

# Add d values we do it by hand because it is easier than parsing the file name...
datbif_meanf <- Map(function(d, df) { data.frame(d = d, df) }, 
                    c(.005, 0.075, 0, 0.15, 0.2, 0.4), 
                    datbif_meanf)

datbif_meanf <- do.call(rbind, datbif_meanf)    

# We drop the second root when d == 0
datbif_meanf <- subset(datbif_meanf, !( d == 0 & root == "root2"))

# Root2 is always the stable one, or root1 if alone
datbif_meanf[ ,'eq_type'] <- ifelse(datbif_meanf[ ,'root'] == "root2", 
                                    "stable", "unstable")
datbif_meanf[datbif_meanf[ ,'d'] == 0,'eq_type'] <- "stable"


# We keep only three values out of the five that highlight what's in the comment 
#   in the pdf doc. 
datbif_meanf <- subset(datbif_meanf, d %in% c(0, 0.15, 0.4))
datbif_meanf[ ,'hetero_rank'] <- rank2(datbif_meanf[ ,'d'])


# Merge both datasets
datbif_musselbed <- rbind.fill(data.frame(type = "Mean-field approx.", datbif_meanf), 
                               data.frame(type = "Spatial simulation", 
                                          eq_type = "stable", # we only have stable eqs in 
                                          datbif_spatial))
                     
datbif_musselbed[ ,'homo'] <- datbif_musselbed[ ,'delta']
datbif_musselbed[ ,'hetero'] <- datbif_musselbed[ ,'d']
datbif_musselbed[ ,'pretty_hetero'] <- with(datbif_musselbed, paste0("d = ", d))
datbif_musselbed[ ,'model'] <- "Musselbed model"
datbif_musselbed[ ,'branch'] <- NA # for compat with other graphs






# Merge all three together
# ------------------

cols <- c('model', 'type', 'homo', 'hetero', 'hetero_rank',
          'pretty_hetero', 'rhop', 'eq_type', 'branch')
datbif_all <- rbind(datbif_grazing[ ,cols],
                    datbif_forestgap[ ,cols], 
                    datbif_musselbed[ ,cols])

# Build a df for the heterogeneous stressor values. Some parameters hand-picked 
#   with love
label_df <- unique(datbif_all[ ,c('model', 'type', 'hetero_rank', 'pretty_hetero')])
label_df[ ,'x'] <- .75
label_df[ ,'y'] <- ifelse(label_df[ ,'type'] == "Mean-field approx.", 
                          4 - .3 * (as.numeric(label_df[ ,'hetero_rank']) - 1),
                          1 - .070 * (as.numeric(label_df[ ,'hetero_rank']) - 1))


ggplot(datbif_all) + 
  geom_line(aes(x = homo, y = rhop, 
                linetype = eq_type, color = hetero_rank,
                group = paste(eq_type, branch, hetero_rank))) + 
  geom_text(aes(x = x, y = y, color = hetero_rank, label = pretty_hetero), 
            hjust = 0, vjust = 1, data = label_df) + 
  facet_grid(type ~ model, switch = "y", scales = 'free_y') + 
  ylab( expression(rho[symbol("+")]) ) + 
  xlab( "Homogeneous stressor" ) + 
  scale_color_manual(values = c('#34AA35','#C9B217','#E8643A')) +
  theme_ifcam + no_legend

ggsave(filename = paste0(output_figure_path, "bifurc_graph_allmodels.pdf"), 
       height = 5*1.5, width = 9*1.5,
       title = 'Bifurcation diagrams for the three models')

