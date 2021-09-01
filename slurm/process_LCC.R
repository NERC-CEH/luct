## ----startup, eval=TRUE, echo=FALSE-------------------------------------------
## set projections and common extent, read in lookup table for reclassification
here::i_am("slurm/process_CORINE.R")
library(spCEH)
library(data.table)
library(here)

# read arguments
# which year in v_times to process
i   <- as.numeric(commandArgs(trailingOnly = TRUE))[1]
print(i)

# # r_25m <- getRasterTemplate(domain = "UK", res = 25, crs = crs_OSGB)

## ----classification, eval=FALSE, echo=FALSE-----------------------------------
# read the LCC_to_LULUCF lookup into a data table
# # dt_u_lcc_lulucf <- fread("../data-raw/LCC/LCC_to_LULUCF.csv")


## ----Reclassify, eval=FALSE, echo=FALSE---------------------------------------
## define the function to rasterize and reclassify LCC data to required classification and extent
reclassify_LCC <- function(year, dt_u_lcc_lulucf, r_25m){

  print(paste0(Sys.time(),": Processing LCC ",year,"..."))

  ## read in LCC data. NB 2015 data has different filename on download

  if(y==2015){

    v_U <- sf::st_read(dsn = paste0("./",year), layer = "cropmap_2015_with_status_170315")

  }else{

    v_U <- sf::st_read(dsn = paste0("./",year), layer = paste0("crops",year,"_GB"))

  }

    st_crs(v_U) <- suppressWarnings(crs(r_25m))

    # standardise the column name
    names(v_U)[grep(names(v_U),pattern = "crop")] <- "crop"

    # merge in the LULUCF codes directly onto the vector data
    v_U_lcc <- dplyr::left_join(v_U, dt_u_lcc_lulucf, by = c("crop"="LCCROP"),copy = F)

    # rasterize the vector data, set 0 to NA and write it
    r_U_lcc <- fasterize(sf = v_U_lcc, raster = r_25m, field = "LULUCF_class_ID")
    writeRaster(r_U_lcc, paste0("./data-raw/LCC/Level1/r_U_lcc_25m_",year,".tif"), overwrite=T)

    print(paste0(Sys.time(),": COMPLETE"))
    return(r_U_lcc)

}



## ----process, eval=TRUE, echo=FALSE-------------------------------------------
v_times <- c(2015, 2016, 2017, 2018, 2019)
v_fnames <- paste0("r_U_lcc_25m_", v_times, ".tif")
v_fnames <- here("data-raw/LCC/Level1", v_fnames)
file.exists(v_fnames)
i <- as.numeric(commandArgs(trailingOnly = TRUE))[1]
year <- v_times[i]

# # r_U_lcc_25m <- reclassify_LCC <- function(year, dt_u_lcc_lulucf, r_25m)

# # fname <- here("data/LCC/Level1", 
  # # paste0("r_U_lcc_25m_", year, ".tif"))
# # writeRaster(r_U_lcc_25m,  fname,  overwrite=T)

## ----agg_disagg, eval=TRUE, echo=FALSE----------------------------------------
# read in re-classified data
r_U_lcc_25m <- raster(v_fnames[i])

# 100 m
r_U_lcc_100m <- aggregate(r_U_lcc_25m, fact = 4, fun = modal)
fname <- here("data/LCC/Level1", 
  paste0("r_U_lcc_100m_", year, ".tif"))
writeRaster(r_U_lcc_100m, fname, overwrite=T)

# 1 km
r_U_lcc_1km  <-  aggregate(r_U_lcc_25m, fact = 1000/25, fun = modal)
fname <- here("data/LCC/Level1", 
  paste0("r_U_lcc_1000m_", year, ".tif"))
writeRaster(r_U_lcc_1km,  fname,  overwrite=T)

# 10 km
r_U_lcc_10km <-  aggregate(r_U_lcc_25m,  fact = 10000/25, fun = modal)
fname <- here("data/LCC/Level1", 
  paste0("r_U_lcc_10000m_", year, ".tif"))
writeRaster(r_U_lcc_10km, fname, overwrite=T)

quit(save = "no")

