
## ----startup, eval=TRUE, echo=FALSE-------------------------------------------
# conditional compilation flags
evalMin    <- TRUE   # always evaluate these chunks, to produce minimla html output
evalSome   <- TRUE  # evaluate if running a subset of chunks
evalAll    <- FALSE  # evaluate if running the whole thing
evalSLURM  <- TRUE  # evaluate if using SLURM
## Load packages & set wd, also install github repository for standardised extents for BELUC model
here::i_am("slurm/process_CORINE.R")
library(spCEH)
library(data.table)
library(here)

# read arguments
# which year in v_times to process
i   <- as.numeric(commandArgs(trailingOnly = TRUE))[1]
print(i)

# set projections and template rasters
## set projections and common extent, read in lookup table for reclassification
# moved to spCEH
#projTM75 <- CRS("+init=epsg:29903") # Irish Grid
r_25m    <- getRasterTemplate(domain = "UK", res = 25, crs = crs_OSGB)
r_NI_25m <- getRasterTemplate(domain = "Northern Ireland", res = 25, crs = crs_OSGB)


## ----classification, eval=evalMin, echo=FALSE---------------------------------
# read the LCM_to_LULUCF lookup into a data table 
fname <- here("data-raw/LCM_WPA", "LCM_to_LULUCF.csv")
dt_u_lcm_lulucf <- fread(fname)
#names(dt_u_lcm_lulucf)
df1 <- dt_u_lcm_lulucf[, 1:4]
knitr::kable(df1)


## ----reclassify_function, eval=evalMin, echo=TRUE-----------------------------
## define the function to reclassify LCM data to required classification and extent
reclassify_LCM <- function(year, dt_u_lcm_lulucf, r_25m){
  
  ## N.B. This process was tried in a cluster over 12 nodes using a UK raster split into 12 parts
  ## but the RAM usage was very high in raster::subs and it wasnt efficient
  print(paste0(Sys.time(),": Processing LCM",year,"..."))
  print(paste0(Sys.time(),":     Combining GB and NI LCMs into UK LCM..."))

  # read in LCM surfaces for GB & NI
  # GB:
  fname <- here(paste0("data-raw/LCM_WPA/",year,"/gb"), 
                paste0("gb",year,"lcm25m.tif"))
  r_U_GB <- raster(fname)
  crs(r_U_GB) <- crs(r_25m)
  # NI:
  fname <- here(paste0("data-raw/LCM_WPA/",year,"/ni"), 
              paste0("ni",year,"lcm25m.tif"))
  r_U_NI_tm75 <- raster(fname)
  crs(r_U_NI_tm75) <- crs_Ire
  
  # NI needs reprojecting to BNG
  r_U_NI <- projectRaster(r_U_NI_tm75, r_NI_25m, method = "ngb")
  
  # merge together into a UK surface and write for posterity & future processing
  r_U_lcm <- merge(r_U_GB, r_U_NI, ext = r_25m)
    fname <- here(paste0("data-raw/LCM_WPA/",year,"/uk"), 
              paste0("uk",year,"lcm25_BNG.tif"))
  writeRaster(r_U_lcm, fname, overwrite=T)
  
  ## substitute values in LCM for new classification (~50 mins) and set 0 to NA
  print(paste0(Sys.time(),":     Reclassifying and writing UK LCM..."))
  
  r_U_lcm <- subs(r_U_lcm, dt_u_lcm_lulucf, by=1, which=3)
  
  r_U_lcm[r_U_lcm == 0] <- NA

  # mask out sea and non-UK land - write new LULUCF raster
  r_U_lcm <- maskByCountry(r_U_lcm, c("England", "Scotland", "Wales", "Northern Ireland"))
  # writeRaster(r_U_lcm, paste0("./data-raw/LCM/Level1/r_U_lcm_25m_",year,".tif"), overwrite=T)

  print(paste0(Sys.time(),": COMPLETE"))
  return(r_U_lcm)
}


## ----process, eval=evalSome, echo=TRUE----------------------------------------
v_times <- c(1990, 2015, 2017, 2018, 2019)
v_fnames <- paste0("r_U_lcm_25m_", v_times, ".tif")
v_fnames <- here("data-raw/LCM_WPA/Level1", v_fnames)
file.exists(v_fnames)
#i =1
year <- v_times[i]


## ----reclassify, eval=evalAll, echo=TRUE--------------------------------------
r_U_lcm_25m <- reclassify_LCM(year = year, dt_u_lcm_lulucf, r_25m)
## write out re-classified UK raster
fname <- here("data/LCM/Level1", 
  paste0("r_U_lcm_25m_", year, ".tif"))
writeRaster(r_U_lcm_25m,  fname,  overwrite=T)

## ----agg_disagg, eval=evalSome, echo=TRUE-------------------------------------
# read in re-classified data
r_U_lcm_25m <- raster(v_fnames[i])

# 100 m
system.time(
  r_U_lcm_100m <- aggregate(r_U_lcm_25m, fact = 4, fun = modal)
)
fname <- here("data/LCM/Level1", 
  paste0("r_U_lcm_100m_", year, ".tif"))
writeRaster(r_U_lcm_100m, fname, overwrite=T)

# 1 km
r_U_lcm_1km  <-  aggregate(r_U_lcm_25m, fact = 1000/25, fun = modal)
fname <- here("data/LCM/Level1", 
  paste0("r_U_lcm_1000m_", year, ".tif"))
writeRaster(r_U_lcm_1km,  fname,  overwrite=T)

# 10 km
r_U_lcm_10km <-  aggregate(r_U_lcm_25m,  fact = 10000/25, fun = modal)
fname <- here("data/LCM/Level1", 
  paste0("r_U_lcm_10000m_", year, ".tif"))
writeRaster(r_U_lcm_10km, fname, overwrite=T)


## ----quit_gracefully, eval=evalSLURM, echo=FALSE------------------------------
# it used to help if scripts exited R without saving when run on LSF job scheduler
# guessing still true with SLURM job scheduler
quit(save = "no")
