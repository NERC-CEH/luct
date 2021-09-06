## ----render_book, eval=FALSE, echo=FALSE, warning=FALSE, message=FALSE--------
## module load jaspy/3.7/r20190627
## #install.packages("bookdown")
## bookdown::render_book("index.Rmd")
## bookdown::preview_chapter("FC.Rmd")
## knitr::purl("FC.Rmd")
## source("slurm/process_FC.R")
#R CMD BATCH --no-restore --no-save slurm/process_FC.R data/FC/log/console.Rout &


## ----startup, eval=TRUE, echo=FALSE, warning=FALSE, message=FALSE-------------
## Load packages 
#source("../R/luc_track.R")
here::i_am("slurm/process_FC.R")
library(targets)
library(tictoc) # capture execution time
library(spCEH)
library(data.table)
library(here)
library(readxl)
library(fasterize)
library(units)
library(rgdal)
library(sf)
library(stars)
# start the execution time clock
tictoc::tic("Computation time (excl. render)")

l_GL_FC  <- tar_read(c_GL_FC)
df_G <- l_GL_FC$df_G
#rasterOptions(tmpdir = "/work/scratch-nopw/bkruijt/raster")


## ----rasterise, eval=TRUE, echo=FALSE-----------------------------------------
# R script to rasterise FC National Forest Estate & Woodland GB data using the R raster package.
# See http://cran.r-project.org/web/packages/raster/index.html
#
# Peter Levy, CEH Edinburgh

r_25    <- getRasterTemplate("UK", res = 25,   crs = crs_OSGB)
r_100   <- getRasterTemplate("UK", res = 100,  crs = crs_OSGB)
r_1000  <- getRasterTemplate("UK", res = 1000, crs = crs_OSGB)
r_10000 <- getRasterTemplate("UK", res = 10000, crs = crs_OSGB)

# beginCluster()
# ogrInfo("../data-raw/FC/National_Forest_Inventory_Woodland_GB_2018-shp")
# ogrInfo("../data-raw/FC/National_Forest_Estate_Subcompartments_England_2019-shp")
# ogrInfo("../data-raw/FC/National_Forest_Estate_Subcompartments_Scotland_2019-shp")
# ogrInfo("../data-raw/FC/NRW_PRODUCTISED_SCDB_LLE")
# ogrInfo("../data-raw/FSNI/FSNI_subcompartments")

# read in NFI GB shapefile
sply_nfiw_gb   <- readOGR(here("data-raw/FC", "National_Forest_Inventory_Woodland_GB_2018-shp"))
dim(sply_nfiw_gb)
levels(sply_nfiw_gb$CATEGORY)
# remove non-woodland polygons
sply_nfiw_gb <- subset(sply_nfiw_gb, CATEGORY == "Woodland")

# NFI does not contain planting dates, so we guess them based on the distribution used in CARBINE
# Distribution will be ok, but spatially will be wrong
# Strictly, should do this only on the ones without a NFES entry, where planting date is known
# Alternatively, NFI might contain all the old forest, planted pre-1900
# so could say uniform sampling 1000-1500 or some other arbitrary date range in the past
# NFI contains 84 % of the area ever planted according to CARBINE
sum(sply_nfiw_gb$Area_ha) / 100 # ha -> km2
sum(df_G$area)

v_sampleYears <- sample(
  x = df_G$time,                          # choosing from this number of years
  size = length(sply_nfiw_gb$CATEGORY), # choosing this number of years
  replace = TRUE, 
  prob = df_G$area)           # with the area planted as the weight
sply_nfiw_gb$PRI_PLYEAR <- v_sampleYears
hist(sply_nfiw_gb$PRI_PLYEAR)
summary(sply_nfiw_gb)

dim(sply_nfiw_gb)
projection(sply_nfiw_gb) <- projOSGB
sf_nfiw_gb <- st_as_sf(sply_nfiw_gb)
# rasterize the polygons
system.time(
  r_U_nfiw_gb_10000 <- fasterize(sf_nfiw_gb, r_10000, field = "PRI_PLYEAR")
)
system.time(
  r_U_nfiw_gb_1000 <- fasterize(sf_nfiw_gb, r_1000, field = "PRI_PLYEAR")
)
system.time(
  r_U_nfiw_gb_100  <- fasterize(sf_nfiw_gb, r_100, field = "PRI_PLYEAR")
)
system.time(
  r_U_nfiw_gb_25   <- fasterize(sf_nfiw_gb, r_25, field = "PRI_PLYEAR")
)

# writeRaster(r_U_nfiw_gb_1000, file = here("data/FC/Level1", "r_U_nfiw_gb_1000.tif"), overwrite = TRUE)
# writeRaster(r_U_nfiw_gb_100,  file = here("data/FC/Level1", "r_U_nfiw_gb_100.tif"), overwrite = TRUE)

sply_nfes_scot <- readOGR(here("data-raw/FC", "National_Forest_Estate_Subcompartments_Scotland_2019-shp"))
sply_nfes_engl <- readOGR(here("data-raw/FC", "National_Forest_Estate_Subcompartments_England_2019-shp"))
sply_nfes_wale <- readOGR(here("data-raw/FC", "NRW_PRODUCTISED_SCDB_LLE"))
sply_nfes_nire <- readOGR(here("data-raw/FSNI", "FSNI_subcompartments"))
projection(sply_nfes_scot) <- projOSGB # proj4 string differs slightly
projection(sply_nfes_engl) <- projOSGB # proj4 string differs slightly
projection(sply_nfes_wale) <- projOSGB # proj4 string differs slightly
# readOGR will read the .prj file if it exists; extract it
sply_nfes_nire <- spTransform(sply_nfes_nire, projOSGB) 

sply_nfes_scot <- sply_nfes_scot[, c("PRI_SPCODE", "PRI_PLYEAR", "PRI_YIELD")]
sply_nfes_engl <- sply_nfes_engl[, c("PRI_SPCODE", "PRI_PLYEAR", "PRI_YIELD")]
sply_nfes_wale <- sply_nfes_wale[, c("PRI_SPCODE", "PRI_PLYEAR", "PRI_YIELD")]
# change names to match
sply_nfes_nire$PRI_SPCODE <- sply_nfes_nire$DESCRIPTIO 
sply_nfes_nire$PRI_PLYEAR <- sply_nfes_nire$YEAR_PLANT 
sply_nfes_nire$PRI_YIELD  <- NA
sply_nfes_nire <- sply_nfes_nire[, c("PRI_SPCODE", "PRI_PLYEAR", "PRI_YIELD")]

# to add the above together:
sply_nfes <- rbind(sply_nfes_scot, sply_nfes_engl, sply_nfes_wale, sply_nfes_nire)
# convert to an sf object
sf_nfes <- st_as_sf(sply_nfes)

# recode PRI_SPCODE from string to unique integers
sf_nfes$PRI_SPCHAR <- sf_nfes$PRI_SPCODE
sf_nfes$PRI_SPCODE <- as.integer(sf_nfes$PRI_SPCODE)
sf_nfes$PRI_PLYEAR <- as.integer(sf_nfes$PRI_PLYEAR)

table(sf_nfes$PRI_PLYEAR)
# recode missing values
sf_nfes$PRI_PLYEAR[sf_nfes$PRI_PLYEAR == 0] <- NA
sf_nfes$PRI_PLYEAR[sf_nfes$PRI_PLYEAR == 9999] <- NA
sf_nfes$PRI_YIELD [sf_nfes$PRI_YIELD == 0]  <- NA
hist(sf_nfes$PRI_PLYEAR)

# sample planting year in proportion
v_sampleYears <- sample(
  x = df_G$time,                           # choosing from these years
  size = sum(is.na(sf_nfes$PRI_PLYEAR)), # choosing this number of years
  replace = TRUE,                       
  prob = df_G$area)            # with the area planted as the weight
hist(v_sampleYears)
sf_nfes$PRI_PLYEAR[is.na(sf_nfes$PRI_PLYEAR)] <- v_sampleYears
hist(sf_nfes$PRI_PLYEAR)

# NFES contains only 30 % of the area ever planted according to CARBINE
# cf. NFI contains   84 % of the area ever planted according to CARBINE
set_units(sum(st_area(sf_nfes)), km^2)
sum(df_G$area)

# rasterize the polygons
system.time(
  r_PlYear_10000m <- fasterize(sf_nfes, r_10000, field = "PRI_PLYEAR")
)
system.time(
  r_PlYear_1000m  <- fasterize(sf_nfes, r_1000, field = "PRI_PLYEAR")
)
system.time(
  r_PlYear_100m   <- fasterize(sf_nfes, r_100, field = "PRI_PLYEAR")
)
system.time(
  r_PlYear_25m    <- fasterize(sf_nfes, r_25, field = "PRI_PLYEAR")
)

r_U_nfiw_gb_1000
compareRaster(r_U_nfiw_gb_1000, r_PlYear_1000m)
cellStats(!is.na(r_U_nfiw_gb_1000), sum)
cellStats(!is.na(r_PlYear_1000m), sum)

# add NFI data where NFES is NA and write to file
r_PlYear_10000m <- cover(r_PlYear_10000m, r_U_nfiw_gb_10000)
r_PlYear_1000m  <- cover(r_PlYear_1000m, r_U_nfiw_gb_1000)
r_PlYear_100m   <- cover(r_PlYear_100m, r_U_nfiw_gb_100)
#r_PlYear_25m   <- cover(r_PlYear_25m, r_U_nfiw_gb_25)

writeRaster(r_PlYear_10000m,  filename = here("data/FC/Level1", "r_PlYear_10000m.tif"), overwrite=TRUE)
writeRaster(r_PlYear_1000m,   filename = here("data/FC/Level1", "r_PlYear_1000m.tif"), overwrite=TRUE)
writeRaster(r_PlYear_100m,    filename = here("data/FC/Level1", "r_PlYear_100m.tif"), overwrite=TRUE)
#writeRaster(r_PlYear_25m,   filename = here("data/FC/Level1", "r_PlYear_25m.tif"), overwrite=TRUE)

# # timing test
# system.time(raster::writeRaster(r_PlYear_100m, filename = here("data/FC/Level1", "r_PlYear_100m.tif"), overwrite=TRUE))
# system.time(terra::writeRaster(r_PlYear_100m, filename = here("data/FC/Level1", "r_PlYear_100m.tif"), overwrite=TRUE))
# system.time(stars::write_stars(st_PlYear_100m, dsn = here("data/FC/Level1", "st_PlYear_100m.tif"), type = "Byte"))

# # convert to stars to test files writing speed
# st_PlYear_1000m <- st_as_stars(r_PlYear_1000m)
# st_PlYear_100m <- st_as_stars(r_PlYear_100m)
# # writing to TIF or netCDF via stars seems x10 quicker; files are x5.5 bigger
# write_stars(st_PlYear_1000m, dsn = here("data/FC/Level1", "st_PlYear_1000m.nc"), type = "Int16")
# write_stars(st_PlYear_100m,  dsn = here("data/FC/Level1", "st_PlYear_100m.nc"), type = "Int16")
# # writing to netCDF via raster seems very slow
# rf <- writeRaster(rNfew_f, filename="rNfew_PRI_PLYEAR_100m.nc",  format="CDF", overwrite=TRUE)

# pdf("r_PlYear_1000m.pdf")
  # plot(st_PlYear_1000m, col=rainbow(25), main="Planting Year, 1000-m resolution")
# dev.off()

# what is the forest area in km2?
cellStats(!is.na(r_U_nfiw_gb_1000), sum)
cellStats(!is.na(r_PlYear_1000m), sum)

#endCluster()

tictoc::toc()

quit(save = "no")
