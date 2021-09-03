## ----start_up, eval=TRUE, echo=FALSE, warning=FALSE, message=FALSE------------
# functions for processing IACS data are country/year dependent and are read in from a separate R script
#source("../R/luc_track.R")
#source("../R/IACS_process_functions.R")
#rasterOptions(tmpdir = "/work/scratch-nopw/bkruijt/raster")

## Load packages 
here::i_am("slurm/process_IACS.R")
library(spCEH)
library(data.table)
library(here)

i <- as.numeric(commandArgs(trailingOnly = TRUE))[1]
v_times <- 2005:2019
v_fnames <- here("data-raw/IACS/Level1", 
  paste0("r_U_iacs_25m_E_0.5_", v_times, ".tif"))
year <- v_times[i]
file.exists(v_fnames)

# set projections and template rasters
r_25m    <- getRasterTemplate(domain = "UK", res = 25, crs = crs_OSGB)

## ----agg_disagg, eval=TRUE, echo=FALSE----------------------------------------
# read in re-classified data
r_U_iacs_25m <- raster(v_fnames[i])

# 100 m
r_U_iacs_100m <- aggregate(r_U_iacs_25m, fact = 100/25, fun = modal)
fname <- here("data/IACS/Level1", 
  paste0("r_U_iacs_100m_", year, ".tif"))
writeRaster(r_U_iacs_100m, fname, overwrite=T)

# 1 km
r_U_iacs_1km  <-  aggregate(r_U_iacs_25m, fact = 1000/25, fun = modal)
fname <- here("data/IACS/Level1", 
  paste0("r_U_iacs_1000m_", year, ".tif"))
writeRaster(r_U_iacs_1km,  fname,  overwrite=T)

# 10 km
r_U_iacs_10km <-  aggregate(r_U_iacs_25m,  fact = 10000/25, fun = modal)
fname <- here("data/IACS/Level1", 
  paste0("r_U_iacs_10000m_", year, ".tif"))
writeRaster(r_U_iacs_10km, fname, overwrite=T)

quit(save = "no")
