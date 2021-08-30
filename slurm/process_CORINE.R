## ----start_up, eval=TRUE, echo=FALSE, warning=FALSE, message=FALSE------------
# For serial run without SLURM
# R CMD BATCH --no-restore --no-save process_CORINE.R console.Rout &
# C:/Progra~1/R/R-4.0.4/bin/R CMD BATCH --no-restore --no-save "--args 1 10000" slurm/process_CORINE.R console.Rout

# read arguments
# which year in v_times to process
i   <- as.numeric(commandArgs(trailingOnly = TRUE))[1]
# resolution of the output
res <- as.numeric(commandArgs(trailingOnly = TRUE))[2]
print(i)
print(res)
#res <- 10000

here::i_am("slurm/process_CORINE.R")
library(spCEH)
library(data.table)
library(here)
# this will load a host of packages used across processing scripts, in one location.
# pkgs <- c("units", "tidyr", "sp", "sf", "raster", "rgeos", "rgdal", "grid", "data.table", "spCEH", "scico", "ggplot2", "ggforce", "plyr", "dplyr", "stars", "BayesianTools", "ggthemes")
# pkgs_available <- sapply(pkgs, require, character.only = TRUE)
# if (!all(pkgs_available)) stop(paste("Not all packages available - check: ", which(!pkgs_available)))
rasterOptions(tmpdir = "/work/scratch-nopw/bkruijt/raster")


## ----Reclassify, eval=TRUE, echo=FALSE----------------------------------------
# Reclassify CORINE LAND COVER
# define the function to reclassify LCm data to required classification and extent
reclassify_CORINE <- function(year, fname_cor, dt_u_corine_lulucf, r_100m_laea, r_100m){
  
  print(paste0(Sys.time(),": Processing CORINE", year,"..."))
  print(paste0(Sys.time(),":     Extracting and reprojecting CORINE data for UK..."))
  
  # read in CORINE data
  r_U_cor_laea <- raster(fname_cor)

  # CORINE needs clipping to reprojected UK extent
  r_U_cor_laea <- crop(r_U_cor_laea, extent(r_100m_laea))
  r_U_cor <- projectRaster(from = r_U_cor_laea, to = r_100m, res = res(r_100m), method = "ngb")
  # fname <- paste0("../data-raw/CORINE/r_Ucor_cor_100m_",year,".tif")
  # writeRaster(r_U_cor, fname, overwrite=T)   
  # }
  
  ## substitute values in CORINE for new classification (~50 mins) and set other values to NA
  print(paste0(Sys.time(),":     Reclassifying and writing UK CORINE..."))
  
  r_U_cor <- subs(r_U_cor, dt_u_corine_lulucf, by=1, which=3)
  
  r_U_cor[r_U_cor == 0] <- NA
  r_U_cor[r_U_cor > 44] <- NA
  
  # mask out sea and non-UK land
  r_U_cor <- maskByCountry(r_U_cor, c("England", "Scotland", "Wales", "Northern Ireland"))
  print(paste0(Sys.time(),": COMPLETE"))
  return(r_U_cor)
}


## ----process, eval=TRUE, echo=FALSE, warning=FALSE, message=FALSE-------------
# set projections and template rasters
projLAEA <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"
projLAEA <- CRS(SRS_string='EPSG:3035')

r_100m <- getRasterTemplate(domain = "UK", res = res, crs = crs_OSGB)
r_100m_laea <- projectRaster(from = r_100m, crs=projLAEA, res=c(res, res))

# read the CORINE_to_LULUCF lookup into a data table
fname <- here("data-raw/CORINE", "CORINE_to_LULUCF.csv")
dt_u_corine_lulucf <- fread(fname)

# Processing CORINE
# For each available CORINE year, run the function.
v_times <- c(1990, 2000, 2006, 2012, 2018)
nYears <- length(v_times)

# get convoluted CORINE raw file name
# make a small lookup table with the years of latest update for each year of data (file structure)
dt_cor_dates <- data.table(clc_year = v_times, update = c(2000, 2006, 2012, 2018, 2018))

v_fname <- paste0("u",dt_cor_dates[clc_year == v_times, update],
  "_clc", v_times, "_v2020_20u1_raster100m/DATA/U",
  dt_cor_dates[clc_year == v_times, update],
  "_CLC",v_times,"_V2020_20u1.tif")
v_fname_cor <- here("data-raw/CORINE", v_fname)
file.exists(v_fname_cor)

# Get job index = i in v_times
#i = 1
year      <- v_times[i]
fname_cor <- v_fname_cor[i]
r_U_cor_100m <- reclassify_CORINE(year = year, fname_cor, dt_u_corine_lulucf, r_100m_laea, r_100m)
## write out re-classified UK raster
fname <- paste0("r_U_cor_100m_", year, ".tif")
fname <- here("data/CORINE/Level1", fname)
writeRaster(r_U_cor_100m, fname, overwrite=T)   

## ----agg_disagg, eval=TRUE, echo=FALSE----------------------------------------

# # read in re-classified data
# r_U_cor_100m <- raster(v_fnames[[1]])

# # 1 km
# r_U_cor_1km  <-  aggregate(r_U_cor_100m, fact = 10, fun = modal)
# fname <- paste0("../data-raw/CORINE/Level1/r_U_cor_1km_", year, ".tif")
# writeRaster(r_U_cor_1km,  fname,  overwrite=T)
# # 10 km
# r_U_cor_10km <-  aggregate(r_U_cor_100m,  fact = 100, fun = modal)
# fname <- paste0("../data-raw/CORINE/Level1/r_U_cor_10km_", year, ".tif")
# writeRaster(r_U_cor_10km, fname, overwrite=T)
# # 25 m
# r_U_cor_25m <- disaggregate(r_U_cor_100m, fact = 4)
# fname <- paste0("../data-raw/CORINE/Level1/r_U_cor_25m_", year, ".tif")
# writeRaster(r_U_cor_25m, fname, overwrite=T)

quit(save = "no")


