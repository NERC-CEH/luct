rm(list = ls(all = TRUE))
here::i_am("./analysis/split_u_by_soiltype.R")
library(targets)
library(here)
library(qs)
library(data.table)
library(stringi)
# library(rasterVis)
# library(ggthemes)
source(here("R/luc_track.R"))
source(here("R/luct.R"))

v_times <- 1950:2020
v_region <- c("en", "sc", "wa", "ni", "uk")
region <- "en"
i_region <- match(region, v_region) # int 1:4
dir_output <- paste0("output/output_", region)
v_soil_type = c("min", "org")
soil_type = "min"
i_soil_type <- match(soil_type, v_soil_type) - 1 # int 0:1

res <- 1000
i_sample <- 1

# read mask dt to land and add region variable
fname <- here("data/mask", paste0("dt_land_", res, "m.qs"))
dt_mask <- qread(fname)
dim(dt_mask)
dt_mask

# read U output file
# fname <- here("output/output_en/dt_u_en_smp1_1000m.qs")
fname <- here(dir_output, paste0("dt_u_", region, "_smp", i_sample, "_", res, "m.qs"))
dt <- qread(fname)
dt
dt <- data.table(dt_mask, dt)

# subset to region and soil_type of interest
dt <- dt[country == i_region & soil_type == i_soil_type]
dim(dt)
dt
#plot(dt$x, dt$y)

# condense into unique land-use vector format
dt_luv <- dt[, .N, by = u_ch]
dt_luv[, area := N * res^2 / 1e6] # res in m^2, area in km^2
dt_luv[, start_time := v_times[1]]
dt_luv[, end_time := v_times[length(v_times)]]
dt_luv <- dt_luv[order(-area)]
dt_luv <- dt_luv[!is.na(u_ch)] # remove NA, only first row with all non-land grid cells
sum(dt_luv[, area])
dim(dt_luv)
# save output
fname <- here(dir_output, paste0("dt_luv_", region, "_", soil_type, "_smp", i_sample, "_", res, "m.qs"))
qsave(dt_luv, file = fname)
