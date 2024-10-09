source("R/luc_track.R")
n_u <- length(names_u)
res <- "10km"

# Read in Level 1 LCM data
v_times <- c(1990, 2015, 2017, 2018, 2019)
v_fnames <- here("data-raw/LCM/Level1",
           paste0("r_U_lcm_", res, "_", v_times, ".tif"))
s_U <- raster::stack(v_fnames)
st_U <- st_as_stars(s_U)
st_U <- st_set_dimensions(st_U, "band", values = v_times, names = "time")
dim(st_U[[1]])
plot(st_U)
st_U_lcm <- st_U

# CORINE
# Read in Level 1 CORINE data
v_times <- c(2000, 2006, 2012, 2018)
v_fnames <- here("data-raw/CORINE/Level1",
           paste0("r_U_cor_", res, "_", v_times, ".tif"))
s_U <- stack(v_fnames)
st_U <- st_as_stars(s_U)
st_U <- st_set_dimensions(st_U, "band", values = v_times, names = "time")
st_U_cor <- st_U

# LCC
# Read in Level 1 Land Cover Crop data
v_times <- c(2015, 2016, 2017, 2018, 2019)
v_fnames <- here("data-raw/LCC/Level1",
           paste0("r_U_lcc_", res, "_", v_times, ".tif"))
s_U <- stack(v_fnames)
st_U <- st_as_stars(s_U)
st_U <- st_set_dimensions(st_U, "band", values = v_times, names = "time")
st_U_lcc <- st_U

# IACS
# Read in Level 1 Land Cover Crop data
v_times <- 2005:2014
v_fnames <- here("data-raw/IACS/Level1",
           paste0("r_U_iacs_", res, "_E_0.5_", v_times, ".tif"))
s_U <- stack(v_fnames)
st_U <- st_as_stars(s_U)
st_U <- st_set_dimensions(st_U, "band", values = v_times, names = "time")
# plot(st_U)
st_U_iac <- st_U
#** need to extend to uk raster from England


# initialise_predictions
# better - declare a star for likelihood with dims x,y,time,u
v_times <- 1990:2019 # that we want to calculate likelihood at
n_t <- length(v_times)
x <- st_get_dimension_values(st_U_lcm, "x")
y <- st_get_dimension_values(st_U_lcm, "y")
time <- v_times
dims = st_dimensions(x = x, y = y, time = time, u = 1:6, cell_midpoints = c(T, T, F, F))
v_dims = dim(dims)
L = array(rep(NA, prod(v_dims)), dim = v_dims)
#** should be a stars proxy object, as it is 781 GB
st_L = st_as_stars(list(L = L), dimensions = dims, proxy = TRUE)
st_dimensions(st_L)
st_L <- st_set_dimensions(st_L, "u", values = names_u)
st_crs(st_L) <- st_crs(st_U_lcm)

# estimate_L
prec_lcm <- 5
prec_lcc <- 1
prec_cor <- 1
prec_iac <- 1
v_prec <- c(prec_lcm, prec_lcc, prec_cor)

# define function to get the layer from a star corresponding to i_time
get_vU_fromStar_byTime <- function(st_U, i_time){
  n_s <- dim(st_U[[1]])[1] * dim(st_U[[1]])[2]
  v_times <- st_get_dimension_values(st_U, "time")

  # find where this year is in stars object
  i_t <- findInterval(i_time, v_times)
  if (i_t > 0) {st_U %>% slice(time, i_t) %>%
    pull() %>% as.vector() -> v_U
  } else {
    v_U <- rep(NA, n_s)
  }
  return(v_U)
}

for (i_t in 1:n_t){
#i_t = 28
  i_time <-  v_times[i_t]
  v_U_lcm <- get_vU_fromStar_byTime(st_U_lcm, i_time)
  v_U_lcc <- get_vU_fromStar_byTime(st_U_lcc, i_time)
  v_U_cor <- get_vU_fromStar_byTime(st_U_cor, i_time)
  #v_U_iac <- get_vU_fromStar_byTime(st_U_iac, i_time)
  #order has to match v_prec i.e. prec_lcm, prec_lcc, prec_cor
  dt <- data.table(v_U_lcm, v_U_lcc, v_U_cor)

  for (i_u in 1:n_u){
  #i_u = 3
    # does each data set observe this land use?
    dt_u <- dt[, lapply(.SD, "==", i_u)]
    # multiply by precision of each data source
    dt_u <- dt_u[, Map("*", .SD, v_prec)]
    #table(dt, useNA = "always")
    #table(dt_u)
    v_L <- rowMeans(dt_u, na.rm = TRUE)
    #summary(v_L)
    #table(v_L, useNA = "always")
    # assign this to the array data from a stars object
    # "[[" retrieves the array from a stars object
    # "[" references the array dimensions (x,y,t,u)
    st_L$L[,, i_t, i_u] <- v_L
    #l_b_L_u[[i_u]][[i_t]] <- setValues(l_b_L_u[[i_u]][[i_t]], v_L)
  } # n_u
}   # n_t

plot(st_L["L",,, c(1,27), 4])
plot(st_L[,,, i_t, i_u])

# stars cannot write out 4-D array
# need to write out 6 layers separately
#** try ncdfgeom or gdalcubes for this
st_L_spl <- split(st_L, "u")
names(st_L_spl) <- paste0("L_u", 1:6)
#lapply(1:6, function(i) write_stars(st_L_spl[names(st_L_spl)[i]], paste0("proxy_test_st_L_u", i, "_", res, ".tif")))

v_fnames <- paste0("st_L_u", 1:n_u, "_", res, ".tif")
for (i in 1:6){
  write_stars(st_L_spl[names(st_L_spl)[i]], v_fnames[i])
}

# DA3plotLloop
v_times <- c(2019)
# v_fnames <- paste0("../data/LikU/st_L_u", 1:n_u, "_", res, ".tif")
v_fnames <- paste0("st_L_u", 1:n_u, "_", res, ".tif")
v_fnames <- paste0("../data/LikU_SPEED/", res, "/s_L_", v_times, ".tif")
v_fnames <- paste0(res, "/st_L_", v_times, ".tif")
st_L <- read_stars(v_fnames)
names(st_L) <- v_times
st_L <- st_set_dimensions(st_L, "band", values = names_u, names = "u")

# for (i in 1:n_u){
  # r <- as(st_L[i,,, c(1)], "Raster")
  # col <- diverge.color(r*100, start.color = "white", end.color = colour_u[i])
  # plot(st_L[i,,, c(1,30)]*100,
    # main = c(paste(names_u[i], "1990"), paste(names_u[i], "2019")),
    # breaks = 0:505, col = col)
  #plot(r, col = col)
# }
```

<!--- } -->

<!--- { DA3plotL -->
```{r DA3plotL1, eval=TRUE, echo=FALSE, warning=FALSE, message=FALSE, fig.cap = paste("Spatial variation in the likelihood $\\mathcal{L}$ of observing land use $u$ in 2019, arbitrarily rescaled for plotting.")}
i <- 1
r <- as(st_L[,,,i], "Raster")
col <- diverge.color(r*100, start.color = "white", end.color = colour_u[i])
plot(st_L[i,,, c(1,30)]*100, main = c(paste(names_u[i], "1990"), paste(names_u[i], "2019")), breaks = 0:505, col = col)
plot(st_L[,,,i]*100, main = c(paste(names_u[i], "2019")), breaks = 0:105, col = col)
```
