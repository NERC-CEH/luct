'%!in%' <- function(x,y)!('%in%'(x,y))

names_u_chess <- c("white", "empty", "black")
names_u  <- c("woods", "crops", "grass", "rough", "urban", "other")
colour_u <- c("lightblue", "green", "pink", "purple", "orange", "midnightblue")

## ----startup, eval=TRUE, echo=FALSE, warning=FALSE, message=FALSE---------------------
# library(plyr)
# library(dplyr)
# library(purrr)
# library(units)
# library(data.table)
# library(ggplot2)

# library(ggforce)
# source("../R/luc_track.R")


#' Function to remove units attribute from matrix or array objects 
#' Units attributes prevent some functions from working
#' 
#' @param x A vector, matrix or array.
#' @return An object of the unchanged class but without units attribute.
#' @export
#' @examples
#' a_B <- remove_units(a_B)
remove_units <- function(x){
    y = as.vector(x)
    dim(y) = dim(x)
    return(y)
}

#' Function to calculate the gross gain in area 
#'  for each land use from a land-use transition matrix
#'
#' @param v_B A land-use transition matrix, or its vector form
#' @param n_u The number of land-use classes
#' @return A vector of the gross gain in area by land use
#' @export
#' @examples
#' A_gain <- getAreaGrossGain_fromBeta(v_B, n_u)
getAreaGrossGain_fromBeta <- function(v_B, n_u = sqrt(length(v_B))){
  m_B <- matrix(v_B, n_u, n_u)
  diag(m_B) <- 0
  A_gain <- colSums(m_B)
  return(A_gain)
}

#' Function to calculate the gross loss in area 
#'  for each land use from a land-use transition matrix
#'
#' @param v_B A land-use transition matrix, or its vector form
#' @return A vector of the gross loss in area by land use
#' @param n_u The number of land-use classes
#' @export
#' @examples
#' A_loss <- getAreaGrossLoss_fromBeta(v_B, n_u)
getAreaGrossLoss_fromBeta <- function(v_B, n_u = sqrt(length(v_B))){
  m_B <- matrix(v_B, n_u, n_u)
  diag(m_B) <- 0
  A_loss <- rowSums(m_B)
  return(A_loss)
}

#' Function to calculate the net change in area 
#'  for each land use from a land-use transition matrix
#'
#' @param v_B A land-use transition matrix, or its vector form
#' @param n_u The number of land-use classes
#' @return A vector of the net change in area by land use
#' @export
#' @examples
#' D_u <- getAreaNetChange_fromBeta(v_B, n_u = 6)
getAreaNetChange_fromBeta <- function(v_B, n_u = sqrt(length(v_B))){
  m_B <- matrix(v_B, n_u, n_u)
  diag(m_B) <- 0
  D_u <- colSums(m_B) - rowSums(m_B)
  return(D_u)  
}

## ---- wrangle_CS

#' Function to wrangle Countryside Survey data 
#'  from text files to R objects
#'
#' @param fpath Filepath to text file
#' @return A BLAG object
#' @export
#' @examples
#' blag <- wrangle_CS(fpath = "../data-raw/CS/UK_LUC_matrices_2018i.csv")
wrangle_CS <- function(fpath = "../data-raw/CS/UK_LUC_matrices_2018i.csv"){
  df_B <- read.csv(fpath)
  # names(df_B)
  # summary(df_B)
  df_B[,6:35] <- df_B[,6:35] * 1000 # raw data are in kha
  # set units to hectares
  df_B[,6:35] <- df_B[,6:35] %>% map2_dfc("ha", ~set_units(.x, .y, mode = "standard"))
  # and convert units to km2
  df_B[,6:35] <- df_B[,6:35] %>% map2_dfc("km^2", ~set_units(.x, .y, mode = "standard"))
  m_B <- df_B[,6:35]
  # dim(m_B)
  # names(m_B)


  ## ----restructure, eval=TRUE, echo=FALSE, warning=FALSE, message=FALSE-----------------
  df_B$year <- df_B$Year
  df_B$Year <- NULL
  #df_B <- subset(df_B, year <= 2020)
  v_times <- unique(df_B$year)
  n_t     <- length(v_times)
  n_u <- length(names_u)
  a_B_ijt <- array(0, c(n_u, n_u, n_t))
  from <- rep(c(1, 4, 3, 2, 5, 6), each = 5)
  to <- c(   4, 3, 2, 5, 6,
          1,    3, 2, 5, 6,
          1, 4,    2, 5, 6,
          1, 4, 3,    5, 6,
          1, 4, 3, 2,    6,
          1, 4, 3, 2, 5   )
  #length(from); length(to)
  n_elem <- length(m_B) # no. elements with entries in matrix; diag is missing
  v_regions <- levels(as.factor(df_B$SubRegion))
  n_regions <- length(v_regions) # no. regions in data frame
  l_a_B <- rep(list(a_B_ijt), n_regions)
  # str(l_a_B)
  n_rows <- dim(df_B)[1]

  for (i_row in 1:n_rows){
  #i_row = 102
    i_time   <- df_B$year[i_row] # 1950 to 2050
    i_t <- match(i_time, v_times) # 1 to 101
    i_region <- df_B$SubRegion[i_row]  # E, NI, S, W
    i_region <- match(i_region, v_regions) # 1:4
    for (i_elem in 1:n_elem) { # col 1:30
    #i_elem = 1
      i_u <- from[i_elem] # from land use
      j_u <- to[i_elem]   # to   land use
      l_a_B[[i_region]][i_u,j_u,i_t] <- m_B[i_row, i_elem]
    }   # i_elem
  }     # itime

  # sum to give the matrices of UK totals
  a_B_uk <- Reduce("+", l_a_B)
  # str(a_B_uk)
  # l_a_B[[1]][,, 71] +
  # l_a_B[[2]][,, 71] +
  # l_a_B[[3]][,, 71] +
  # l_a_B[[4]][,, 71]
  # a_B_uk[,, 71]


  ## ----CSplotBtext, eval=TRUE, echo=FALSE, warning=FALSE, message=FALSE, fig.cap = "Example $\\mathbf{B}$ matrix showing the areas changing land use between 1990 and 1991 in km^2^."----
  # B matrix
  a_B <- a_B_uk
  #st_B <- st_as_stars(a_B)
  #plot(st_B)
  a_B <- remove_units(a_B)
  #par(mfrow=c(2,2))
  #apply(a_B[,,2:5], MARGIN=3, function(x) {corrplot(x, is.corr = FALSE, insig = "p-value")})
  #plot(st_B, text_values = TRUE)

  dt_B <- as.data.table(a_B[,,]) # remove first time - no difference to calculate, so all zeroes
  names(dt_B) <- c("u_from", "u_to", "time", "area")
  # summary(dt_B)
  #dt_B$area <- remove_units(dt_B$area)
  dt_B$year <- v_times[dt_B$time]
  dt_B$u_from <- as.factor(dt_B$u_from)
  dt_B$u_to   <- as.factor(dt_B$u_to)
  levels(dt_B$u_from) <- names_u
  levels(dt_B$u_to)   <- names_u
  #dt_B <- subset(dt_B, year > 1970 & year < 2020)
  # dt_B <- subset(dt_B, year == 1970 | 
    # year == 1980 | year == 1990 | year == 2000)

  ## ----CSplotG, eval=TRUE, echo=FALSE, warning=FALSE, message=FALSE, fig.cap = "Time series of implied area gains $\\mathbf{G}$ to each land use, from CS data. We assumed that the rates of change were constant during the period between surveys."----
  # calc gross gain from B; MARGIN = 3 because that is the time dimension
  dt_G <- as.data.table(t(apply(a_B[,,], MARGIN = 3, FUN = getAreaGrossGain_fromBeta, n_u = 6)))
  names(dt_G) <- names_u
  dt_G <- data.table(time = as.numeric(rownames(dt_G)), dt_G)
  dt_G <- melt(dt_G, id=c("time"), variable.name = "u", value.name = "area")
  dt_G$year <- v_times[dt_G$time]

  ## ----CSplotL, eval=TRUE, echo=FALSE, warning=FALSE, message=FALSE, fig.cap = "Time series of implied losses of area $\\mathbf{L}$ from each land use, from CS data. We assumed that the rates of change were constant during the period between surveys."----
  # calc gross loss from B; MARGIN = 3 because that is the time dimension
  dt_L <- as.data.table(t(apply(a_B[,,], MARGIN = 3, FUN = getAreaGrossLoss_fromBeta, n_u = 6)))
  names(dt_L) <- names_u
  dt_L <- data.table(time = as.numeric(rownames(dt_L)), dt_L)
  dt_L <- melt(dt_L, id=c("time"), variable.name = "u", value.name = "area")
  dt_L$year <- v_times[dt_L$time]

  ## ----CSplotD, eval=TRUE, echo=FALSE, warning=FALSE, message=FALSE, fig.cap = "Time series of implied net change in area $\\Delta \\mathbf{A}$ of each land use, from CS data. We assumed that the rates of change were constant during the period between surveys."----
  # dt_D_cs, ## this is D not A**
  # calc net change from B; MARGIN = 3 because that is the time dimension
  dt_D <- as.data.table(t(apply(a_B[,,], MARGIN = 3, FUN = getAreaNetChange_fromBeta, n_u = 6)))
  names(dt_D) <- names_u
  dt_D <- data.table(time = as.numeric(rownames(dt_D)), dt_D)
  dt_D <- melt(dt_D, id=c("time"), variable.name = "u", value.name = "area")
  dt_D$year <- v_times[dt_D$time]

  ## ----saving, eval=TRUE, echo=FALSE, warning=FALSE, message=FALSE----------------------
  dt_D$data_source <- "CS"
  dt_B$data_source <- "CS"
  dt_G$data_source <- "CS"
  dt_L$data_source <- "CS"

  dt_D$time <- dt_D$year; dt_D$year <- NULL 
  dt_B$time <- dt_B$year; dt_B$year <- NULL 
  dt_G$time <- dt_G$year; dt_G$year <- NULL 
  dt_L$time <- dt_L$year; dt_L$year <- NULL

  dt_D$area <- set_units(dt_D$area, km^2)
  dt_B$area <- set_units(dt_B$area, km^2)
  dt_G$area <- set_units(dt_G$area, km^2)
  dt_L$area <- set_units(dt_L$area, km^2)

  dt_D$area <- set_units(dt_D$area, m^2)
  dt_B$area <- set_units(dt_B$area, m^2)
  dt_G$area <- set_units(dt_G$area, m^2)
  dt_L$area <- set_units(dt_L$area, m^2) 

  dt_D$area <- drop_units(dt_D$area)
  dt_B$area <- drop_units(dt_B$area)
  dt_G$area <- drop_units(dt_G$area)
  dt_L$area <- drop_units(dt_L$area) 

  dt_D_cs <- dt_D
  dt_B_cs <- dt_B
  dt_G_cs <- dt_G
  dt_L_cs <- dt_L 
  return(list(dt_D = dt_D, dt_B = dt_B, dt_G = dt_G, dt_L = dt_L))
}


## ---- wrangle_FC

#' Function to wrangle FC time series data 
#'  from text files to R objects
#'
#' @param v_fpath Vector of filepaths to data files
#' @return A list of data frames for affn and defn i.e. gross gains and losses
#' @export
#' @examples
#' fpath1 = "./data-raw/FC/timeSeries/forest_planting_byYear_UK.csv"
#' fpath2 = "./data-raw/FC/timeSeries/Deforestation_Areas_for_CEH_1990-2019.xlsx"
#' x <- wrangle_FC(c(fpath1, fpath2))
wrangle_FC <- function(v_fpath = 
  c("./data-raw/FC/timeSeries/forest_planting_byYear_UK.csv",
    "./data-raw/FC/timeSeries/Deforestation_Areas_for_CEH_1990-2019.xlsx")){

  # read data on afforestation
  df_affn <- read.csv(v_fpath[1])
  names(df_affn)
  names(df_affn) <- c("time", "area")
  df_affn$area <- set_units(df_affn$area, km^2)
  df_affn$area <- set_units(df_affn$area, m^2)
  summary(df_affn)
  # plot(df_affn$time, cumsum(df_affn$area))
  df_affn$u <- "woods"
  df_affn$data_source <- "FC"
  dt_G <- as.data.table(df_affn)
  dt_G$area <- drop_units(dt_G$area)

  # read data on deforestation
  df_defn <- read_excel(v_fpath[2],
    sheet  = "Def Time series 1990-2019i", skip = 3)
  df_defn <- with(df_defn[1:30,], data.frame(time = as.numeric(Year), area = Total))
  df_defn$area <- set_units(df_defn$area, ha)
  df_defn$area <- set_units(df_defn$area, m^2)
  df_defn$u <- "woods"
  df_defn$data_source <- "FC"
  dt_L <- as.data.table(df_defn)
  dt_L$area <- drop_units(dt_L$area)

  # initialise a copy for net change
  dt_D <- dt_L
  dt <- merge(dt_L, dt_G, all.x = TRUE, by = "time")
  dt_D$area <- dt$area.y - dt$area.x
  dt_D$area <- set_units(dt_D$area, m^2)
  dt_D$area <- drop_units(dt_D$area)
  
  # initialise a copy for absolute area
  dt_A  <- dt_D
  # assume initial area of forest in 1990 worked out previously - check this
  initArea <- (36312.07 * 1e6) - sum(dt_D$area)
  dt_A$area <- initArea + cumsum(dt_D$area)
  
  return(list(dt_A = dt_A, dt_D = dt_D, dt_G = dt_G, dt_L = dt_L))
}


## ---- run_corine_job

#' Function to run Corine processing job 
#'  from raw tif files to R objects
#'
#' @param fname_job File path to SLURM job file for CORINE processing
#' @return A job object
#' @export
#' @examples
#' fname_job = "./slurm/process_CORINE.job"
#' x <- run_corine_job(fname_job)
run_corine_job <- function(fname_job = "./slurm/process_CORINE.job"){
  cmd <- paste0("sbatch ", fname_job)
  # submit the jobs to the SLURM queue
  err <- system(cmd)
  # and return the years and paths of the output files
  # these need to match slurm/process_CORINE.R - no programmed check they are consistent
  v_times <- c(2000, 2006, 2012, 2018)
  return(list(
    v_times = v_times,
    v_fnames = paste0("./data/CORINE/Level1/r_U_cor_100m_", v_times, ".tif")
  ))
}


## ---- run_iacs_job

#' Function to run IACS processing job 
#'  from raw tif files to R objects
#'
#' @param fname_job File path to SLURM job file for IACS processing
#' @return A job object
#' @export
#' @examples
#' fname_job = "./slurm/process_IACS.job"
#' x <- run_iacs_job(fname_job)
run_iacs_job <- function(fname_job = "./slurm/process_IACS.job"){
  cmd <- paste0("sbatch ", fname_job)
  # submit the jobs to the SLURM queue
  err <- system(cmd)
  # and return the years and paths of the output files
  # these need to match slurm/process_IACS.R - no programmed check they are consistent
  v_times <- 2005:2019
  return(list(
    v_times = v_times,
    v_fnames = paste0("./data/IACS/Level1/r_U_iacs_100m_", v_times, ".tif")
  ))
}


## ---- run_FC_job

#' Function to run FC processing job 
#'  from raw shp files to R objects
#'
#' @return A log file path
#' @export
#' @examples
#' fname_log <- run_FC_job()
run_FC_job <- function(fname_code = "slurm/process_FC.R"){
  # initialise an empty job object
  cmd <- paste0("R CMD BATCH --no-restore --no-save slurm/", fname_code, " data/FC/log/console.Rout &")
  # submit the jobs and get the time to identify the output files from this batch
  err <- system(cmd)
  # and return the path of the first output file
  return(fname_out = "data/FC/Level1/r_PlYear_100m.tif")  
}


## ---- run_lcc_job

#' Function to run LCC processing job 
#'  from raw tif files to R objects
#'
#' @param fname_job File path to SLURM job file for LCC processing
#' @return A job object
#' @export
#' @examples
#' fname_job = "./slurm/process_LCC.job"
#' x <- run_lcc_job(fname_job)
run_lcc_job <- function(fname_job = "./slurm/process_LCC.job"){
  cmd <- paste0("sbatch ", fname_job)
  # submit the jobs to the SLURM queue
  err <- system(cmd)
  # and return the years and paths of the output files
  # these need to match slurm/process_LCC.R - no programmed check they are consistent
  v_times <- 2015:2019
  return(list(
    v_times = v_times,
    v_fnames = paste0("./data/LCC/Level1/r_U_lcc_100m_", v_times, ".tif")
  ))
}


## ---- run_lcm_job

#' Function to run LCM processing job 
#'  from raw tif files to R objects
#'
#' @param fname_job File path to SLURM job file for LCM processing
#' @return A job object
#' @export
#' @examples
#' fname_job = "./slurm/process_LCM.job"
#' x <- run_lcm_job(fname_job)
run_lcm_job <- function(fname_job = "./slurm/process_LCM.job"){
  # no reprocessing needed now - done by Sam outwith targets
  # redo later within targets
  # cmd <- paste0("sbatch ", fname_job)
  # # submit the jobs to the SLURM queue
  # err <- system(cmd)
  # and return the years and paths of the output files
  # these need to match slurm/process_LCM.R - no programmed check they are consistent
  v_times <- c(1990, 1994, 1998, 2002, 2006, 2010, 2015, 2017, 2018, 2019)
  return(list(
    v_times = v_times,
    v_fnames = paste0("./data/LCM/Level1/r_U_lcm_1000m_", v_times, ".tif")
  ))
}

## ---- apply_mask

#' Function to mask Level1 raster files 
#'  with a consistent land mask
#'
#' @param level1 List with vector of times and file path for Level 1 files
#' @param res Numeric Resolution of raster files in metres
#' @return The same level1 object; side effect is that the files have been re-written
#' @export
#' @examples
#' fname_job = "./slurm/process_LCM.job"
#' x <- apply_mask(level1, res = 100)
apply_mask <- function(level1, res = 100){
  level1$v_fnames <- str_replace(level1$v_fnames, "1000", as.character(res))
  
  n_files <- length(level1$v_fnames)
  n_values_raw    <- vector(mode = "integer", length = n_files)
  n_values_masked <- vector(mode = "integer", length = n_files)
  fname <- paste0("./data-raw/mask/r_land_", res, "m.tif")
  r_mask <- raster(fname)
  # plot(r_mask)
  # plot(r)
  n <- cellStats(!is.na(r_mask), sum)
  
  for (i in 1:n_files){
  #i = 1
    fname <- level1$v_fnames[i] # paste0("./data/LCM/Level1/r_U_lcm_1000m_", v_times, ".tif")
    r <- raster(fname)
    n_values_raw[i] <- cellStats(!is.na(r), sum)
    r <- mask(r, r_mask)
    n_values_masked[i] <- cellStats(!is.na(r), sum)
    writeRaster(r, fname, overwrite = TRUE)
  }
  # plot(n_values_raw / n * 100)
  # points(n_values_masked / n * 100, col = "red")
  
  return(level1)
}
## ---- run_crome_job

#' Function to run CROME processing job 
#'  from raw tif files to R objects
#'
#' @param fname_job File path to SLURM job file for CROME processing
#' @return A job object
#' @export
#' @examples
#' fname_job = "./slurm/process_CROME.job"
#' x <- run_crome_job(fname_job)
run_crome_job <- function(fname_job = "./slurm/process_CROME.job"){
  # no reprocessing needed now - done by Sam outwith targets
  # no job called, just returns the Level1 file paths
  # redo later within targets
  # cmd <- paste0("sbatch ", fname_job)
  # # submit the jobs to the SLURM queue
  # err <- system(cmd)
  # and return the years and paths of the output files
  # these need to match slurm/process_CROME.R - no programmed check they are consistent
  v_times <- c(2016, 2017, 2018, 2019)
  return(list(
    v_times = v_times,
    v_fnames = paste0("./data/CROME/Level1/r_U_crome_100m_", v_times, ".tif")
  ))
}


## ---- combine_blags

#' Function to combine observations
#'  in BLAG objects to produce a single data table
#'
#' @param l_blags List of blag objects to combine
#' @return A blag object
#' @export
#' @examples
#' l_blags = list(c_blag_CS, c_blag_corine, c_blag_iacs, c_blag_lcc, c_blag_lcm)
#' x <- combine_blags(l_blags)
combine_blags <- function(l_blags = list(blag_cor, blag_iacs, blag_lcc, blag_lcm)){
  # B matrix
  dt_B <- rbindlist(lapply(l_blags, '[[', "dt_B"), use.names=TRUE, fill=TRUE)
  
  # re-ordering by (u_to, u_from) (not u_from, u_to) allows 
  # the default vector to matrix conversion to work (where byrow = FALSE)
  dt_B <- arrange(dt_B, time, data_source, u_to, u_from)

  # G time series
  dt_G <- rbindlist(lapply(l_blags, '[[', "dt_G"), use.names=TRUE, fill=TRUE)
  
  # L time series
  dt_L <- rbindlist(lapply(l_blags, '[[', "dt_L"), use.names=TRUE, fill=TRUE)
  
  # D time series
  dt_D <- rbindlist(lapply(l_blags, '[[', "dt_D"), use.names=TRUE, fill=TRUE)

  # A time series
  dt_A <- rbindlist(lapply(l_blags, '[[', "dt_A"), use.names=TRUE, fill=TRUE)
    
  return(list(dt_B = dt_B, dt_G = dt_G, dt_L = dt_L, dt_A = dt_A, dt_D = dt_D))
}


## ---- convert_units

#' Function to convert units of area
#'  in BLAG objects in all the data tables
#'
#' @param blag A blag object
#' @return A blag object
#' @export
#' @examples
#' x <- convert_units(blag, old_unit = "m^2", new_unit = "km^2")
convert_units <- function(blag, old_unit = "m^2", new_unit = "km^2"){
  # set the mode of set_units() to read character variables
  units::units_options(set_units_mode = "standard")

  # B matrix
  blag$dt_B$area <-  set_units(blag$dt_B$area, old_unit)
  blag$dt_B$area <-  set_units(blag$dt_B$area, new_unit)
  blag$dt_B$area <- drop_units(blag$dt_B$area)
  
  # re-ordering by (u_to, u_from) (not u_from, u_to) allows 
  # the default vector to matrix conversion to work (where byrow = FALSE)
  blag$dt_B <- arrange(blag$dt_B, time, data_source, u_to, u_from)

  # G time series
  blag$dt_G$area <-  set_units(blag$dt_G$area, old_unit)
  blag$dt_G$area <-  set_units(blag$dt_G$area, new_unit)
  blag$dt_G$area <- drop_units(blag$dt_G$area)
  
  # L time series
  blag$dt_L$area <-  set_units(blag$dt_L$area, old_unit)
  blag$dt_L$area <-  set_units(blag$dt_L$area, new_unit)
  blag$dt_L$area <- drop_units(blag$dt_L$area)
  
  # D time series
  blag$dt_D$area <-  set_units(blag$dt_D$area, old_unit)
  blag$dt_D$area <-  set_units(blag$dt_D$area, new_unit)
  blag$dt_D$area <- drop_units(blag$dt_D$area)

  # A time series - may not always be present
  if (length(blag$dt_A) > 0){
    blag$dt_A$area <-  set_units(blag$dt_A$area, old_unit)
    blag$dt_A$area <-  set_units(blag$dt_A$area, new_unit)
    blag$dt_A$area <- drop_units(blag$dt_A$area)  
  }

  # set the mode of set_units() back to default
  units::units_options(set_units_mode = "symbols")
  
  return(blag)
}


## ---- set_exclusions

#' Function to exclude observations
#'  in BLAG objects for specific data source x land use combinations
#'
#' @param obs A blag object
#' @param data_sources_toInclude A vector of data sources to include
#' @return A blag object with data sources excluded
#' @export
#' @examples
#' obs_km2 <- tar_read(c_obs_km2)
#' unique(obs$dt_A$data_source[obs$dt_A$u == "other"])
#' obs_exc <- set_exclusions(obs_km2, data_sources_toInclude = 
#'   c("AgCensus", "CS", "FC", "MODIS")) # WP-A subset
#' obs_exc <- set_exclusions(obs_km2)
#' unique(obs$dt_A$data_source[obs$dt_A$u == "other"])
##* WIP obs should be replaced with blag
set_exclusions <- function(
  obs,
  data_sources_toInclude = c("AgCensus", "CORINE", "CROME", "CS", "FC",
                           "IACS",     "LCC",    "LCM",   "MODIS")
  ){
  # subset to only the data sets we want to include
  obs$dt_B <- obs$dt_B[data_source %in% data_sources_toInclude]
  obs$dt_G <- obs$dt_G[data_source %in% data_sources_toInclude]
  obs$dt_L <- obs$dt_L[data_source %in% data_sources_toInclude]
  obs$dt_D <- obs$dt_D[data_source %in% data_sources_toInclude]

  # zeroes are actually missing values
  # should remove earlier
  # obs$dt_A$area[obs$dt_A$area == 0]    <- NA
  obs$dt_D$area[obs$dt_D$area == 0]    <- NA
  obs$dt_D$area[is.nan(obs$dt_D$area)] <- NA
  obs$dt_G$area[obs$dt_G$area == 0]    <- NA
  obs$dt_G$area[is.nan(obs$dt_G$area)] <- NA
  obs$dt_L$area[obs$dt_L$area == 0]    <- NA
  obs$dt_L$area[is.nan(obs$dt_L$area)] <- NA
  obs$dt_B$area[obs$dt_B$area == 0] <- NA

  unique(obs$dt_D$data_source); unique(obs$dt_B$data_source)

  # obs$dt_A$useData <- TRUE
  obs$dt_G$useData <- TRUE
  obs$dt_L$useData <- TRUE
  obs$dt_D$useData <- TRUE

  # exclude these data sources for woods
  data_sources_toExclude <- c("AgCensus", "IACS", "LCC", "CROME")
  # obs$dt_A$useData[obs$dt_A$u == "woods" & obs$dt_A$data_source %in% data_sources_toExclude] <- FALSE
  obs$dt_G$useData[obs$dt_G$u == "woods" & obs$dt_G$data_source %in% data_sources_toExclude] <- FALSE
  obs$dt_L$useData[obs$dt_L$u == "woods" & obs$dt_L$data_source %in% data_sources_toExclude] <- FALSE
  obs$dt_D$useData[obs$dt_D$u == "woods" & obs$dt_D$data_source %in% data_sources_toExclude] <- FALSE

  # exclude these data sources for grass
  data_sources_toExclude <- c("IACS")
  # obs$dt_A$useData[obs$dt_A$u == "grass" & obs$dt_A$data_source %in% data_sources_toExclude] <- FALSE

  # exclude these data sources for rough
  data_sources_toExclude <- c("IACS", "LCC")
  # obs$dt_A$useData[obs$dt_A$u == "rough" & obs$dt_A$data_source %in% data_sources_toExclude] <- FALSE
    
  # exclude these data sources for urban
  # obs$dt_A$useData[obs$dt_A$u == "urban" & obs$dt_A$data_source %in% data_sources_toExclude] <- FALSE
  obs$dt_G$useData[obs$dt_G$u == "urban" & obs$dt_G$data_source %in% data_sources_toExclude] <- FALSE
  obs$dt_L$useData[obs$dt_L$u == "urban" & obs$dt_L$data_source %in% data_sources_toExclude] <- FALSE
    
  # exclude these data sources for other
  # obs$dt_A$useData[obs$dt_A$u == "other" & obs$dt_A$data_source %in% data_sources_toExclude] <- FALSE
  obs$dt_G$useData[obs$dt_G$u == "other" & obs$dt_G$data_source %in% data_sources_toExclude] <- FALSE
  obs$dt_L$useData[obs$dt_L$u == "other" & obs$dt_L$data_source %in% data_sources_toExclude] <- FALSE

  obs$dt_D <- subset(obs$dt_D, useData)
  # obs$dt_A <- subset(obs$dt_A, useData)
  obs$dt_G <- subset(obs$dt_G, useData)
  obs$dt_L <- subset(obs$dt_L, useData)

  return(obs)
}


## ---- get_loglik

#' Define a function to calculate the log-likelihood (loglik) 
#' for observations and a given $B$ matrix. 
#'
#' @param v_B Beta matrix as a vector
#' @return The RMSE from comparison of predictions versus observations
#' @return dt_B_obs Data table containing observations of Beta matrix terms
#' @return dt_G_obs Data table containing observations of gross gains
#' @return dt_L_obs Data table containing observations of gross losses
#' @return dt_D_obs Data table containing observations of net area changes
#' @return n_u Numeric The number of land uses
#' @return use_rmse = FALSE
#' @export
#' @examples
#' obs <- tar_read(c_obs)
#' rmse   <- get_loglik(v_B, use_rmse = TRUE)
#' rmse   <- get_loglik(v_B, obs$dt_B, obs$dt_G, obs$dt_L, obs$dt_D, use_rmse = TRUE)
#' loglik <- get_loglik(v_B, obs$dt_B, obs$dt_G, obs$dt_L, obs$dt_D)
get_loglik <- function(v_B, 
    #dt_B_obs, dt_G_obs, dt_L_obs, dt_D_obs, 
    n_u = sqrt(length(v_B)), use_rmse = FALSE){
  m_B_pred  <- matrix(v_B, n_u, n_u)
  dt_B_pred <- data.table(u_from = rep(names_u, n_u), u_to = rep(names_u, each = n_u), pred  = as.vector(m_B_pred))
  dt_G_pred <- data.table(u = names_u, pred  = getAreaGrossGain_fromBeta(m_B_pred, n_u))
  dt_L_pred <- data.table(u = names_u, pred  = getAreaGrossLoss_fromBeta(m_B_pred, n_u))
  dt_D_pred <- data.table(u = names_u, pred = getAreaNetChange_fromBeta(m_B_pred, n_u))
  
  dt_B <- merge(dt_B_obs, dt_B_pred, by = c("u_from", "u_to"))
  dt_G <- merge(dt_G_obs, dt_G_pred, by = "u")
  dt_L <- merge(dt_L_obs, dt_L_pred, by = "u")
  dt_D <- merge(dt_D_obs, dt_D_pred, by = "u")

  if (use_rmse){
    resid_D <- dt_D[, area - pred]
    resid_G <- dt_G[, area - pred]
    resid_L <- dt_L[, area - pred]
    resid_B <- dt_B[, area - pred]

    RMSE_D <- sqrt(mean(resid_D^2, na.rm = TRUE)) # * wt_D
    RMSE_G <- sqrt(mean(resid_G^2, na.rm = TRUE)) # * wt_G 
    RMSE_L <- sqrt(mean(resid_L^2, na.rm = TRUE)) # * wt_L 
    RMSE_B <- sqrt(mean(resid_B^2, na.rm = TRUE)) # * wt_B 
    
    v_RMSE <- c(RMSE_D, RMSE_B, RMSE_G, RMSE_L)
    # if no data, these will be NaN, which need to be NA
    v_RMSE[is.nan(v_RMSE)] <- NA
    return(sum(v_RMSE, na.rm = TRUE))
  } else {
    # use constant CV for sigma in obs; could use mult_sd <- c(0.05, 0.1, 0.1, 0.4, 0.3, 0.2)
    # v_llik_B <- dnorm(dt_B$area, mean = dt_B$pred, sd = 0.1*abs(dt_B$area), log = T)
    # v_llik_G <- dnorm(dt_G$area, mean = dt_G$pred, sd = 0.1*abs(dt_G$area), log = T)
    # v_llik_L <- dnorm(dt_L$area, mean = dt_L$pred, sd = 0.1*abs(dt_L$area), log = T)
    # v_llik_D <- dnorm(dt_D$area, mean = dt_D$pred, sd = 0.1*abs(dt_D$area), log = T)

    # crude treatment of Fn rate ~ 1/3 positive of Fp rate
    v_llik_B <- dnorm(dt_B$area * (1- dt_B$Fp * 0.66), mean = dt_B$pred, sd = dt_B$sigma*abs(dt_B$area), log = T)
    v_llik_G <- dnorm(dt_G$area * (1- dt_G$Fp * 0.66), mean = dt_G$pred, sd = dt_G$sigma*abs(dt_G$area), log = T)
    v_llik_L <- dnorm(dt_L$area * (1- dt_L$Fp * 0.66), mean = dt_L$pred, sd = dt_L$sigma*abs(dt_L$area), log = T)
    v_llik_D <- dnorm(dt_D$area * (1- dt_D$Fp * 0.66), mean = dt_D$pred, sd = dt_D$sigma_D, log = T)
    loglik <- sum(v_llik_D, v_llik_G, v_llik_L, v_llik_B, na.rm = TRUE)
    return(loglik)
  }
}


## ---- get_pred_ls

#' Define a function to calculate the least-squares predictions 
#' for the $B$ matrix, and associated G, L, and D data tables. 
#'
#' @param v_B Beta matrix as a vector
#' @return A blag object for the least-squares predictions
#' @export
#' @examples
#' obs <- tar_read(c_obs)
#' pred_ls <- get_pred_ls(obs, start_year = 2018, end_year = 2019)
get_pred_ls <- function(
  obs,
  start_year = 2015,
  end_year = 2019){

  v_times <- start_year:end_year
  n_t <- length(v_times)

  # this duplicates reordering in combine_blags
  # but the additional check is probably a good idea
  # as the order determines vector to matrix conversion ordering
  obs$dt_B <- dplyr::arrange(obs$dt_B, time, data_source, u_to, u_from)

  # initialise array for predicted B parameters
  # use LCM as initial values
  v_B <- obs$dt_B[time == 2019 & data_source == "LCM"]$area
  v_B[is.na(v_B)] <- 0
  m_B <- matrix(v_B, 
    nrow = n_u, ncol = n_u, dimnames = list(names_u, names_u))
  a_B_pred <- array(m_B, c(n_u, n_u, n_t))

  # generate initial values
  # using zeroes (no influence on results)
  v_B_ini <- a_B_pred[ , , 1] * 0
  #diag(v_B_ini) <- diag(a_B_pred[ , , 1])
  # using obs + random
  #v_B_ini <- a_B_obs[ , , i_t] + round(runif(9, 0, 20))
  #v_B_ini <- a_B_pred[,,1]
  # initial test
  # get_rmse(v_B_ini)

  # Loop over each year, finding the $B$ matrix which minimises RMSE with the observed data.
  # backwards version
  #for (i_t in (n_t):2){
  # forwards version
  start_time <- Sys.time()
  for (i_t in 2:n_t){
  #i_t = 2
    i_time <- v_times[i_t]
    # need to super-assign these to the global environment
    # because that is where get_loglik is defined, so the local
    # values are not available to it (scoping is lexical not dynamic).
    dt_B_obs <<- obs$dt_B[time == i_time]
    dt_G_obs <<- obs$dt_G[time == i_time]
    dt_L_obs <<- obs$dt_L[time == i_time]
    dt_A_obs <<- obs$dt_A[time == i_time]
    dt_D_obs <<- obs$dt_D[time == i_time]

    # minimise the RMSE
    fit <- optim(v_B_ini, get_loglik, 
      # dt_B_obs = dt_B_obs, dt_G_obs = dt_G_obs, 
      # dt_L_obs = dt_L_obs, dt_D_obs = dt_D_obs,
      use_rmse = TRUE,      
      method = "L-BFGS-B", 
      lower = 0, control = list(factr = 1e9, maxit = 1000, trace = 2))

    # # alternatively, maximise the log-likelihood instead
    # # use_rmse = FALSE (default) and set fnscale = -1 to maximise
    # # seems slower
    # fit <- optim(fit$par, get_loglik, 
      # dt_B_obs = dt_B_obs, dt_G_obs = dt_G_obs, 
      # dt_L_obs = dt_L_obs, dt_D_obs = dt_D_obs,  
      # method = "L-BFGS-B", 
      # lower = 0, control = list(fnscale = -1, factr = 1e9, maxit = 1000, trace = 2))

   # fit <- optim(fit$par, getRMSE, method = "L-BFGS-B", 
     # lower = 0, control = list(factr = 1e3, maxit = 1000, trace = 2))

    # compare B
    v_B_ini
    round(fit$par, 0)
    # a_B_obs[ , , i_t]
    # store B matrix for this time step in array
    a_B_pred[ , , i_t] <- matrix(fit$par, n_u, n_u)
  }

  print(paste0("time elapsed (mins): ", round(Sys.time() - start_time, 3)))

  # derive the Losses And Gains from the Beta matrix 
  pred <- getBLAG(a_B_pred, names_u = names_u, name_data_source = "pred_ls",
    v_times = v_times)  

  # not needed I think
  # # Beta matrix
  # pred$dt_B <- subset(pred$dt_B, time > v_times[1])
  # # net change
  # pred$dt_D <- subset(pred$dt_D, time > v_times[1])
  # # gross Gain
  # pred$dt_G <- subset(pred$dt_G, time > v_times[1])
  # # gross Loss
  # pred$dt_L <- subset(pred$dt_L, time > v_times[1])

  return(pred)
}


## ---- get_post_mcmc_serial

#' Define a function to calculate the least-squares predictions 
#' for the $B$ matrix, and associated G, L, and D data tables. 
#'
#' @param v_B Beta matrix as a vector
#' @return A blag object for the least-squares predictions
#' @export
#' @examples
#' obs     <- tar_read(c_obs)
#' pred_ls <- tar_read(c_pred_ls)
#' l_post <- get_post_mcmc_serial(obs, pred_ls, start_year = 2018, end_year = 2019, n_iter = 1000)
get_post_mcmc_serial <- function(
  obs,
  pred_ls,
  start_year = 2015,
  end_year = 2019,
  n_iter = 10,
  n_chains = 3,
  thin = round(max(1, n_iter/1000))
  ){

  v_times <- start_year:end_year
  n_t <- length(v_times)
  # allocate an empty list to store an mcmc output object for each time step
  l_mcmcOut <- vector("list", n_t)
  
  # this duplicates reordering in combine_blags
  # but the additional check is probably a good idea
  # as the order determines vector to matrix conversion ordering
  obs$dt_B <- dplyr::arrange(obs$dt_B, time, data_source, u_to, u_from)

  ## start BayesianTools MCMC code

  # Prior: use CS from 2019
  # v_B_prior <- dt_B_obs[time == 2019 & data_source == "CS", area]
  # v_B_prior[is.na(v_B_prior)] <- 0
  # v_sd_prior <- 0.1*v_B_prior + 0.1
  ## probably over-complicated
  # prior <- createTruncatedNormalPrior(mean = v_B_prior, sd = v_sd_prior, 
    # lower = rep(0, length(v_B_prior)), upper = 3*max(v_B_prior)+10)
  # or just uniform ** or half-normal
  prior <- createUniformPrior(
    lower = rep(    0, n_u^2), 
    upper = rep(10000, n_u^2))
  setUp <- createBayesianSetup(get_loglik, 
    prior = prior, parallel = FALSE)

  # Initial values for chains 
  # get LS predictions as a starting point 
  v_B_ini <- pred_ls$dt_B[time == 2019, area]
  m_starter <- matrix(rep(v_B_ini, n_chains), nrow = n_chains, byrow = TRUE)
  m_starter[2,] <- rep(0, n_u^2)         # all zeroes
  m_starter[3,] <- runif(n_u^2, 0, 1000) # random

  # initialise array for predicted B parameters
  m_B <- matrix(0, 
    nrow = n_u, ncol = n_u, dimnames = list(names_u, names_u))
  a_B_pred <- array(m_B, c(n_u, n_u, n_t))

  # Loop over each year, finding the $B$ matrix which minimises RMSE with the observed data.
  # backwards version
  #for (i_t in (n_t):2){
  # forwards version
  start_time <- Sys.time()
  for (i_t in 2:n_t){
  #i_t = 2
    i_time <- v_times[i_t]
    # need to super-assign these to the global environment
    # because that is where get_loglik is defined, so the local
    # values are not available to it (scoping is lexical not dynamic).
    dt_B_obs <<- obs$dt_B[time == i_time]
    dt_G_obs <<- obs$dt_G[time == i_time]
    dt_L_obs <<- obs$dt_L[time == i_time]
    dt_A_obs <<- obs$dt_A[time == i_time]
    dt_D_obs <<- obs$dt_D[time == i_time]

    # start BayesianTools MCMC code
    settings <- list(iterations = n_iter,  thin = thin,  nrChains = n_chains, message = TRUE)
    system.time(out <- runMCMC(bayesianSetup = setUp, sampler = "DREAMzs", settings = settings))
    l_mcmcOut[[i_t]] <- out
  }
  print(paste0("time elapsed (mins): ", round(Sys.time() - start_time, 3)))

  return(l_mcmcOut)
}


## ---- get_post_mcmc_parallel

#' Define a function to calculate the least-squares predictions 
#' for the $B$ matrix, and associated G, L, and D data tables. 
#'
#' @param v_B Beta matrix as a vector
#' @return A blag object for the least-squares predictions
#' @export
#' @examples
#' obs     <- tar_read(c_obs)
#' pred_ls <- tar_read(c_pred_ls)
#' get_post_mcmc_parallel <- get_post_mcmc(obs, pred_ls, start_year = 2016, end_year = 2019, n_iter = 1000)
get_post_mcmc_parallel <- function(
  obs,
  pred_ls,
  start_year = 2015,
  end_year = 2019,
  n_iter = 10,
  # we want three processors to each run one chain, with the 3 internal chains 
  n_chains = 1, # on the same core (DREAMz uses 3 internal chains by default)
  n_cores  = 3, # number of cores to use
  thin = round(max(1, n_iter/1000))
  ){

  v_times <- start_year:end_year
  n_t <- length(v_times)
  # allocate an empty list to store an mcmc output object for each time step
  l_mcmcOut <- vector("list", n_t)
  
  # this duplicates reordering in combine_blags
  # but the additional check is probably a good idea
  # as the order determines vector to matrix conversion ordering
  obs$dt_B <- dplyr::arrange(obs$dt_B, time, data_source, u_to, u_from)

  ## start BayesianTools MCMC code

  # Prior: use CS from 2019
  # v_B_prior <- dt_B_obs[time == 2019 & data_source == "CS", area]
  # v_B_prior[is.na(v_B_prior)] <- 0
  # v_sd_prior <- 0.1*v_B_prior + 0.1
  ## probably over-complicated
  # prior <- createTruncatedNormalPrior(mean = v_B_prior, sd = v_sd_prior, 
    # lower = rep(0, length(v_B_prior)), upper = 3*max(v_B_prior)+10)
  # or just uniform ** or half-normal
  prior <- createUniformPrior(
    lower = rep(    0, n_u^2), 
    upper = rep(10000, n_u^2))
  setUp <- createBayesianSetup(get_loglik, 
    prior = prior, parallel = FALSE)

  # Initial values for chains 
  # get LS predictions as a starting point 
  v_B_ini <- pred_ls$dt_B[time == 2019, area]
  m_starter <- matrix(rep(v_B_ini, n_cores), nrow = n_cores, byrow = TRUE)
  m_starter[2,] <- rep(0, n_u^2)         # all zeroes
  m_starter[3,] <- runif(n_u^2, 0, 1000) # random


  # initialise array for predicted B parameters
  m_B <- matrix(0, 
    nrow = n_u, ncol = n_u, dimnames = list(names_u, names_u))
  a_B_pred <- array(m_B, c(n_u, n_u, n_t))

  # Loop over each year, finding the $B$ matrix which minimises RMSE with the observed data.
  # backwards version
  #for (i_t in (n_t):2){
  # forwards version
  start_time <- Sys.time()
  for (i_t in 2:n_t){
  #i_t = 2
    i_time <- v_times[i_t]
    dt_B_obs <- obs$dt_B[time == i_time]
    dt_G_obs <- obs$dt_G[time == i_time]
    dt_L_obs <- obs$dt_L[time == i_time]
    dt_A_obs <- obs$dt_A[time == i_time]
    dt_D_obs <- obs$dt_D[time == i_time]

    # start BayesianTools MCMC code
    settings <- list(iterations = n_iter,  thin = thin,  nrChains = n_chains, message = TRUE)
    system.time(out <- runMCMC(bayesianSetup = setUp, sampler = "DREAMzs", settings = settings))
    l_mcmcOut[[i_t]] <- out
  }
  print(paste0("time elapsed (mins): ", round(Sys.time() - start_time, 3)))

  return(l_mcmcOut)
}


## ---- run_mcmc_job

#' Function to run MCMC processing job 
#'  for Beta matrix
#'
#' @param fname_job File path to SLURM job file for LCM processing
#' @return A job object
#' @export
#' @examples
#' fname_job = "./slurm/process_LCM.job"
#' x <- run_mcmc_beta_job(fname_job, obs)
run_mcmc_beta_job <- function(
    fname_job = "./slurm/run_mcmc_beta.job", 
    obs
  ){
  cmd <- paste0("sbatch ", fname_job)
  # submit the jobs to the SLURM queue
  err <- system(cmd)
  # and return the years and paths of the output files
  # these should match slurm/run_mcmc_beta.R, but not essential - only ones that are tracked
  v_times <- 1991:2019
  v_fnames <- fs::path_rel(here("output", paste0("mcmcB_map", v_times, ".qs")))
  return(v_fnames)
}


## ---- get_nrmse

#' Function to RMSE normalised by range 
#'  for values with mean near zero
#'
#' @param df A data frame containing the variables
#' @param v A character string for the name of the test variable
#' @param v_ref A character string for the name of the reference variable
#' @return Numeric The root-mean-square error
#' @export
#' @examples
#' nrmse <- get_nrmse(df = df, v = "IACS", v_ref = "Ref")
get_nrmse <- function(df = df, v, v_ref = "Ref"){
  v     <- df[[v]]
  v_ref <- df[[v_ref]]
  resid <- v - v_ref
  rmse <- sqrt(mean(resid^2, na.rm = TRUE))
  nrmse <- rmse / 
    (max(v_ref, na.rm = TRUE) - min(v_ref, na.rm = TRUE))
  # if no data, these will be NaN, which need to be NA
  nrmse[is.nan(nrmse)] <- NA
  return(nrmse)
}


## ---- get_r2

#' Function to run MCMC processing job 
#'  for Beta matrix
#'
#' @param df A data frame containing the variables
#' @param v A character string for the name of the test variable
#' @param v_ref A character string for the name of the reference variable
#' @return Numeric The r2 correlation coefficient
#' @export
#' @examples
#' r2 <- get_r2(df = df, v = "CORINE", v_ref = "Ref")
get_r2 <- function(df = df, v, v_ref = "Ref"){
  v     <- df[[v]]
  v_ref <- df[[v_ref]]
  r2 <- cor(v, v_ref, use = "pairwise.complete.obs")
  r2
  # if no data, these will be NaN, which need to be NA
  r2[is.nan(r2)] <- NA
  return(r2)
}


## ---- get_uncert_scaling

#' Function to run MCMC processing job 
#'  for Beta matrix
#'
#' @param obs A blag object containing the observations
#' @param v_names_sources A character vector for the names of the data sources
#' @param cv_AgCensus Numeric Coefficient of variation assumed for AgCensus
#' @return df A data frame containing the scaling variables for uncertainty 
#' @export
#' @examples
#' df_uncert <- get_uncert_scaling(obs_exc,   v_names_sources = 
#'   c("AgCensus", "MODIS", "CS", "FC", "LCM", "CORINE", "LCC", "IACS", "CROME"),
#'   v_interval_length_sources = 
#'   c( 1,          1,       8,    1,    3,     6,        1,     1,      1),
#'   v_start_year_source = 
#'   c( 1750,       2001,    1950, 1900, 1990,  1990,     2015,  2004,   2016),
#'   cv_AgCensus = 0.1
#'   )
#' df_uncert
get_uncert_scaling <- function(obs, 
  v_names_sources = 
  c("AgCensus", "MODIS", "CS", "FC", "LCM", "CORINE", "LCC", "IACS", "CROME"),
  v_interval_length_sources = 
  c( 1,          1,       8,    1,    3,     6,        1,     1,      1),
  v_start_year_source = 
  c( 1750,       2001,    1950, 1900, 1990,  2000,     2015,  2004,   2016),
  cv_AgCensus = 0.1
  ){
  
  dt_D <- obs$dt_D
  dt_D$country <- NULL
  df <- pivot_wider(dt_D, names_from = data_source, values_from = area)
  df <- subset(df, time >= 1990 & time <= 2020)
  df$Ref <- df$AgCensus
  df$Ref[df$u == "woods"] <- df$FC[df$u == "woods"]
  df$Ref[df$u == "urban"] <- df$LCM[df$u == "urban"]

  v_nrmse <- sapply(v_names_sources, get_nrmse, df = df, v_ref = "Ref")
  v_r2   <- sapply(v_names_sources, get_r2,   df = df, v_ref = "Ref")

  df <- data.frame(intvl_lth = v_interval_length_sources, 
    start_year_source = v_start_year_source, 
    NRMSE = v_nrmse, r2 = v_r2, 
    # reduce NRMSE proportional to r2, so that abs and prop measures contribute to sigma weighting
    sigma = v_nrmse * abs(1 - v_r2))

  df <- df[order(df$sigma),]
  # AgCensus and FC form the reference, so are rows 1:2 when ordered
  # guess sigma for these as half the lowest value, which will be row 3 
  # very arbitrary assumption, to be improved upon
  df["AgCensus",]$sigma <- df[3,]$sigma * 0.5
  df["FC",]$sigma       <- df[3,]$sigma * 0.5

  # all values are relative to AGCensus, assuming it has 10 % uncertainty
  df$sigma <- df$sigma / df["AgCensus",]$sigma * cv_AgCensus
  
  # add a term to account for non-annual data
  # each year is a sample from the popn of years included
  # like estimating sd from se with n years in sample
  df$sigma <- df$sigma * sqrt(df$intvl_lth)


  # sigma for D has to be in absolute terms, not relative
  # find the mean magnitude of D by data source
  dt_D_mean <- obs$dt_D[, mean(abs(area), na.rm = TRUE), by = data_source][, .(data_source = data_source, D_mean = V1)]

  ##### one-off calc - can  be moved to a meta notebook
  ##### start
  # D values are average to zero, so need absolute uncert, not relative
  dt_G_mean <- obs$dt_G[, mean(abs(area), na.rm = TRUE), by = data_source][, .(data_source = data_source, G_mean = V1)]
  dt_L_mean <- obs$dt_L[, mean(abs(area), na.rm = TRUE), by = data_source][, .(data_source = data_source, L_mean = V1)]
  dt_GL <- merge(dt_G_mean, dt_L_mean) 
  dt <- merge(dt_D_mean, dt_GL, all.x = TRUE, by = "data_source") 
  # gross fluxes are ~3 times larger than net on average
  GL_to_D <- mean(mean(dt[, G_mean / D_mean], na.rm = TRUE),
    mean(dt[, L_mean / D_mean], na.rm = TRUE))
  ##### end
  ##### one-off calc - can  be moved to a meta notebook
  
  df <- merge(df, dt_D_mean, all.x = TRUE, by.x = "row.names", by.y = "data_source")
  # and calc abs value of sigma; multipier because G & L 3 times larger than D
  df$sigma_D <- sqrt((df$sigma * df$D_mean * GL_to_D)^2 * 2)

  # add values for false positive and neg rates:
  df_fpn <- readRDS("./data/df_fpn.rds")
  df <- merge(df, df_fpn, all.x = TRUE, by.x = "Row.names", by.y = "row.names")
  names(df)[names(df) == "Row.names"] <- "data_source"
  df$Fp[is.na(df$Fp)] <- 0
  df$Fn[is.na(df$Fn)] <- 0
  df$sigma[is.na(df$sigma)] <- median(df$sigma, na.rm = TRUE)

  df <- df[order(df$sigma),]
  return(df)
}


## ---- add_uncert

#' Function to add uncertainties to observations 
#'
#' @param obs A blag object containing the observations
#' @param c_df_uncert A data frame containing the scaling variables for uncertainty 
#' @return obs A blag object containing the observations with uncertainties
#' @export
#' @examples
#' obs_unc <- add_uncert(obs_filled, df_uncert)
add_uncert <- function(obs, df_uncert){
  
  #df_uncert$data_source <- rownames(df_uncert)
  df_uncert <- df_uncert[, c("data_source", "sigma", "sigma_D", "start_year_source", "Fp", "Fn")]
  obs$dt_B <- merge(obs$dt_B, df_uncert, all.x = TRUE, by = "data_source") 
  obs$dt_G <- merge(obs$dt_G, df_uncert, all.x = TRUE, by = "data_source") 
  obs$dt_L <- merge(obs$dt_L, df_uncert, all.x = TRUE, by = "data_source") 
  #obs$dt_A <- merge(obs$dt_A, df_uncert, all.x = TRUE, by = "data_source") 
  obs$dt_D <- merge(obs$dt_D, df_uncert, all.x = TRUE, by = "data_source") 
  
  # scale sigma to account for data interpolated prior to start date for data source
  # if time is before start of data, sigma increase with time^2; otherwise multiplier = 1
  obs$dt_B$sigma_time_term <- pmax(obs$dt_B$start_year_source - obs$dt_B$time + 1, 1)^2
  obs$dt_G$sigma_time_term <- pmax(obs$dt_G$start_year_source - obs$dt_G$time + 1, 1)^2
  obs$dt_L$sigma_time_term <- pmax(obs$dt_L$start_year_source - obs$dt_L$time + 1, 1)^2
  obs$dt_D$sigma_time_term <- pmax(obs$dt_D$start_year_source - obs$dt_D$time + 1, 1)^2
                         
  obs$dt_B$sigma <- obs$dt_B$sigma   * obs$dt_B$sigma_time_term
  obs$dt_G$sigma <- obs$dt_G$sigma   * obs$dt_G$sigma_time_term
  obs$dt_L$sigma <- obs$dt_L$sigma   * obs$dt_L$sigma_time_term
  # note option to rename sigma_D to sigma here, so we can use same terms afterwards
  # but possibly more confusing, as they are different units (abs vs relative)
  obs$dt_D$sigma_D <- obs$dt_D$sigma_D * obs$dt_D$sigma_time_term
    
  obs$dt_B <- obs$dt_B[!is.na(area)]
  obs$dt_G <- obs$dt_G[!is.na(area)]
  obs$dt_L <- obs$dt_L[!is.na(area)]
  obs$dt_D <- obs$dt_D[!is.na(area)]
    
  return(obs)
}

## ----- correct_blag

#' Function to correct B, L, A & G in BLAG objects based on false positives added from add_uncert
#' @param obs A blag object that has already had the function add_uncert applied to it (providing the Fp column)
#' @param df_Fp dataframe returned from get_uncert including false positive values
#' @return obs A blag object containing the updated observations based on the false positive and negative rates
#' @export
#' @examples
#' obs_corr <- correct_blag(data_source_in = "FC", blag = c_obs_unc)
#' obs_corr_uncert <- add_uncert(obs_corr, df_uncert)
correct_blag <- function(data_source_in, blag){

  # subset to a single data source
  #dt_B <- blag$dt_B[data_source == data_source_in]
  dt_B <- subset(blag$dt_B, data_source == data_source_in)

  # only correct if B matrix exists - not all data sources
  if (nrow(dt_B) > 0){
    # apply the correction
    dt_B$area <- dt_B$area * (1 - dt_B$Fp)
    
    v_times <- unique(dt_B$time)
    n_t     <- length(v_times)
    n_u <- length(names_u)
    a_B <- array(0, c(n_u, n_u, n_t))
    
    for (i in 1:length(v_times)){
    #i = 3
      v_B <- dt_B[time == v_times[i], area]
      a_B[, , i] <- matrix(v_B, n_u, n_u)
    }

    # calc gross gain from B; MARGIN = 3 because that is the time dimension
    dt_G <- as.data.table(t(apply(a_B[,,], MARGIN = 3, FUN = getAreaGrossGain_fromBeta, n_u = 6)))
    names(dt_G) <- names_u
    dt_G <- data.table(time = as.numeric(rownames(dt_G)), dt_G)
    dt_G <- melt(dt_G, id=c("time"), variable.name = "u", value.name = "area")
    dt_G$year <- v_times[dt_G$time]

    # calc gross loss from B; MARGIN = 3 because that is the time dimension
    dt_L <- as.data.table(t(apply(a_B[,,], MARGIN = 3, FUN = getAreaGrossLoss_fromBeta, n_u = 6)))
    names(dt_L) <- names_u
    dt_L <- data.table(time = as.numeric(rownames(dt_L)), dt_L)
    dt_L <- melt(dt_L, id=c("time"), variable.name = "u", value.name = "area")
    dt_L$year <- v_times[dt_L$time]

    # dt_D_cs, ## this is D not A**
    # calc net change from B; MARGIN = 3 because that is the time dimension
    dt_D <- as.data.table(t(apply(a_B[,,], MARGIN = 3, FUN = getAreaNetChange_fromBeta, n_u = 6)))
    names(dt_D) <- names_u
    dt_D <- data.table(time = as.numeric(rownames(dt_D)), dt_D)
    dt_D <- melt(dt_D, id=c("time"), variable.name = "u", value.name = "area")
    dt_D$year <- v_times[dt_D$time]

    dt_D$data_source <- data_source_in
    dt_G$data_source <- data_source_in
    dt_L$data_source <- data_source_in

    dt_D$time <- dt_D$year; dt_D$year <- NULL 
    dt_G$time <- dt_G$year; dt_G$year <- NULL 
    dt_L$time <- dt_L$year; dt_L$year <- NULL
    
    dt_B <- as.data.table(dt_B)
    dt_G <- as.data.table(dt_G)
    dt_L <- as.data.table(dt_L)
    dt_D <- as.data.table(dt_D)
    return(list(dt_D = dt_D, dt_B = dt_B, dt_G = dt_G, dt_L = dt_L))
  } else {
    return(blag)
  }
}

#' @examples
#' v_data_source = c("AgCensus", "MODIS", "CS", "FC", "LCM", "CORINE", "LCC", "IACS", "CROME")
#' v_data_source_to_correct = c( "CORINE", "CROME", "CS", "IACS", "LCC", "LCM" )
#' v_data_source <- c("IACS", "LCC", "CROME")
#' obs_corr <- get_correct_blag(v_data_source, blag = c_obs_unc)
 get_correct_blag <- function(v_data_source, blag){
  l_blags <- lapply(v_data_source, FUN = correct_blag, blag = blag)
  blag <- combine_blags(l_blags)
  return(blag)
}


## ----- interpolate_blag

#' Function to interpolate B, L, A & G in BLAG objects to remove missing values
#' @param obs A blag object that has already had the function add_uncert applied to it (providing the Fp column)
#' @param start_year Start year for filling NA values 
#' @param end_year   End year for filling NA values 
#' @return obs A blag object containing the interpolated data
#' @export
#' @examples
#' blag <- obs
#' obs_filled <- interpolate_blag(obs_exc, start_year = 1950, end_year = 2020)
interpolate_blag <- function(blag, start_year = 1950, end_year = 2020){
  # define generic functions
  # G, L, & D grouped by u
  # B grouped by u_from and u_to
  
  fill_na_B <- function(dt){
    # subset to time period
    dt <- dt[time >= start_year & time <= end_year]
    # remove previously interpolated CS data
    dt <- dt[!(data_source == "CS" & time %% 10 != "0")]
    # convert time to a factor with all levels 
    dt <- dt[, time := factor(time, levels = start_year:end_year)]
    # add all the missing factor combinations
    dt <- as.data.table(complete(dt, data_source, u_from, u_to, time))
    dt <- dt[, time := as.integer(as.character(time))]
    # interpolate missing values
    dt <- dt[, area:= na.approx(area, na.rm = FALSE), by = .(data_source, u_from, u_to)]
    # fill in trailing missing values
    dt <- dt[, area:= na.locf(area, na.rm = FALSE, fromLast = TRUE), by = .(data_source, u_from, u_to)]
    # fill in leading missing values
    dt <- dt[, area:= na.locf(area, na.rm = FALSE, fromLast = FALSE), by = .(data_source, u_from, u_to)]
    return(dt)
  }

  # make a generic function
  fill_na_GLD <- function(dt){
    # subset to time period
    dt <- dt[time >= start_year & time <= end_year]
    # remove previously interpolated CS data
    dt <- dt[!(data_source == "CS" & time %% 10 != "0")]
    # convert time to a factor with all levels 
    dt <- dt[, time := factor(time, levels = start_year:end_year)]
    # add all the missing factor combinations
    dt <- as.data.table(complete(dt, data_source, u, time))
    dt <- dt[, time := as.integer(as.character(time))]
    # interpolate missing values
    dt <- dt[, area:= na.approx(area, na.rm = FALSE), by = .(data_source, u)]
    # fill in trailing missing values
    dt <- dt[, area:= na.locf(area, na.rm = FALSE, fromLast = TRUE), by = .(data_source, u)]
    # fill in leading missing values
    dt <- dt[, area:= na.locf(area, na.rm = FALSE, fromLast = FALSE), by = .(data_source, u)]
    return(dt)
  }

  blag$dt_B <- fill_na_B(blag$dt_B)
  blag$dt_G <- fill_na_GLD(blag$dt_G)
  blag$dt_L <- fill_na_GLD(blag$dt_L)
  blag$dt_D <- fill_na_GLD(blag$dt_D)
  
  # dim(blag$dt_B)
  # summary(blag$dt_B)
  
  return(blag)
}


## ----- get_post_plots

#' Function to interpolate B, L, A & G in BLAG objects to remove missing values
#' @param obs A blag object that has already had the function add_uncert applied to it (providing the Fp column)
#' @param start_year Start year for filling NA values 
#' @param end_year   End year for filling NA values 
#' @return obs A blag object containing the interpolated data
#' @export
#' @examples
#' blag <- obs
#' obs_filled <- get_post_plots(obs_exc, start_year = 1950, end_year = 2020)

get_post_plots <- function(
  start_time = 1950,
  end_time   = 2020,
  v_times = start_time:end_time,
  dir_output = "output",
  v_mcmc_fname_Bmap = paste0(dir_output, "/mcmcB_", v_times, ".qs"),
  fig_start_time = 1990,
  fig_end_time   = 2020,
  obs_unc = obs_unc,
  obs_exc = obs_exc,
  v_data_source = unique(obs_exc$dt_D$data_source),
  blag_lcm = blag_lcm,
  rethin = 10,
  start  = 9000
  ){
  
  obs_exc$dt_B <- obs_exc$dt_B[data_source %in% v_data_source]
  obs_exc$dt_G <- obs_exc$dt_G[data_source %in% v_data_source]
  obs_exc$dt_L <- obs_exc$dt_L[data_source %in% v_data_source]
  obs_exc$dt_D <- obs_exc$dt_D[data_source %in% v_data_source]
  
  n_t <- length(v_times)
  fname <- here(dir_output, "mcmcB_2019.qs")
  out <- qread(fname)
  # Combine the chains
  out <- createMcmcSamplerList(out)
  #str(out, max.level = 2)
  # extract a sample matrix just to get dimensions
  sample <- getSample(out, start = start, thin = rethin)
  n_s <- dim(sample)[1] # num samples
  names_q <- c("2.5%", "50%", "97.5%")
  n_q <- length(names_q)  # num quantiles = c(0.025, 0.5, 0.975)

  # declare arrays to store predictions
  a_sample    <- array(data = NA, dim = c(n_s, n_u^2, n_t))
  a_B_post    <- array(data = NA, dim = c(n_u, n_u, n_t, n_s))
  a_B_map    <- array(data = NA, dim = c(n_u, n_u, n_t))
  a_A_post    <- array(data = 0, dim = c(n_u, n_s, n_t),
    dimnames = list(names_u, 1:n_s, v_times))
  a_A_post_map <- array(data = 0, dim = c(n_u, n_t),
    dimnames = list(names_u, v_times))
  # set last element to observed area:
  a_A_post[,,n_t] <- matrix(rep(as.numeric(blag_lcm$dt_A_wide[time == 2019, 1:6]), n_s), ncol=n_s, nrow=n_u)
  a_A_post_map[,n_t] <- as.numeric(blag_lcm$dt_A_wide[time == 2019, 1:6])
  a_A_post_q  <- array(data = NA, dim = c(n_q, n_u, n_t),
    dimnames = list(names_q, names_u, v_times))
  a_D_post     <- array(data = NA, dim = c(n_u, n_s, n_t))
  a_D_post_q   <- array(data = NA, dim = c(n_q, n_u, n_t),
    dimnames = list(names_q, names_u, v_times))
  a_D_post_map <- array(data = NA, dim = c(n_u, n_t),
    dimnames = list(names_u, v_times))
  # posterior estimates of Area gross gain and loss
  a_G_post     <- array(data = NA, dim = c(n_u, n_s, n_t))
  a_L_post     <- array(data = NA, dim = c(n_u, n_s, n_t))
  # quantiles for gains and losses
  a_G_post_q   <- array(data = NA, dim = c(n_q, n_u, n_t),
    dimnames = list(names_q, names_u, v_times))
  a_L_post_q   <- array(data = NA, dim = c(n_q, n_u, n_t),
    dimnames = list(names_q, names_u, v_times))
  a_G_post_map <- array(data = NA, dim = c(n_u, n_t),
    dimnames = list(names_u, v_times))
  a_L_post_map <- array(data = NA, dim = c(n_u, n_t),
    dimnames = list(names_u, v_times))
  # quantiles for B parameter
  names_u_matrix <- paste(rep(names_u, n_u), "_to_", rep(names_u, each = n_u), sep="")
  # needed to set facet_wrap order as we want, otherwise does it byrow 
  names_u_matrix_byrow <- paste(rep(names_u, each = n_u), "_to_", rep(names_u, n_u), sep="")

  a_B_post_q   <- array(data = NA, dim = c(n_q, n_u^2, n_t),
    dimnames = list(names_q, names_u_matrix, v_times))
  a_B_post_map <- a_B_post_q[2,,]

  for (t in 2:n_t){
  #t = 29
    fname <- here(dir_output, paste0("mcmcB_", v_times[t], ".qs"))
    out <- qread(fname)
    
    # Combine the chains
    out <- createMcmcSamplerList(out)

    a_sample[,,t] <- getSample(out, start = start, thin = rethin)
    for (s in 1:n_s){ # must be a vectorised way of doing this - re-order dimensions?
      a_B_post[,,t,s] <- matrix(a_sample[s,,t], n_u, n_u)
    }
    a_D_post_q[,,t] <- getPredictiveIntervals(a_sample[,,t], getAreaNetChange_fromBeta, quantiles = c(0.025, 0.5, 0.975))$posteriorPredictiveCredibleInterval
    # get the Maxiumum APosteriori parameter vector
    v_B_map <- MAP(out)$parametersMAP
    # and save MAP prediction values to array
    a_B_map[,,t] <- matrix(v_B_map, n_u, n_u)
    a_B_post_map[,t] <- v_B_map
    # and MAP predictions
    a_D_post_map[,t] <- getAreaNetChange_fromBeta(v_B_map)
    
    a_D_post[,,t] <- apply(a_sample[,,t], FUN = getAreaNetChange_fromBeta, MARGIN = 1)
    a_A_post[,,t-1] <-  a_A_post[,,t] - a_D_post[,,t]
    a_A_post_map[,t-1] <-  a_A_post_map[,t] - a_D_post_map[,t]

    a_G_post_q[,,t] <- getPredictiveIntervals(a_sample[,,t], getAreaGrossGain_fromBeta, quantiles = c(0.025, 0.5, 0.975))$posteriorPredictiveCredibleInterval
    a_L_post_q[,,t] <- getPredictiveIntervals(a_sample[,,t], getAreaGrossLoss_fromBeta, quantiles = c(0.025, 0.5, 0.975))$posteriorPredictiveCredibleInterval
    # and MAP predictions
    a_G_post_map[,t] <- getAreaGrossGain_fromBeta(v_B_map)
    a_L_post_map[,t] <- getAreaGrossLoss_fromBeta(v_B_map)
    # and quantiles for B parameters
    a_B_post_q[,,t] <- getCredibleIntervals(a_sample[,,t], quantiles = c(0.025, 0.5, 0.975))
  }

  # processing MCMC chain output
  # quantiles of area prediction from net change
  df_D_post_q025 <- data.frame(v_times, t(a_D_post_q[1,,]))
  df_D_post_q50  <- data.frame(v_times, t(a_D_post_q[2,,]))
  df_D_post_q975 <- data.frame(v_times, t(a_D_post_q[3,,]))
  df_D_post_map  <- data.frame(v_times, t(a_D_post_map))
                        # melt(dt_G, id=c("time", "dtime"), variable.name = "u", value.name = "area")
  df_D_post_q025_long <- melt(df_D_post_q025, id=c("v_times"), variable.name = "u")
  df_D_post_q50_long  <- melt(df_D_post_q50 , id=c("v_times"), variable.name = "u")
  df_D_post_q975_long <- melt(df_D_post_q975, id=c("v_times"), variable.name = "u")
  df_D_post_map_long  <- melt(df_D_post_map,  id=c("v_times"), variable.name = "u")

  names(df_D_post_q025_long)[3] <- "area_q025"
  names(df_D_post_q50_long)[3]  <- "area_q50"
  names(df_D_post_q975_long)[3] <- "area_q975"
  names(df_D_post_map_long)[3]  <- "area_MAP"

  df_D_post_long <- merge(df_D_post_q025_long, df_D_post_q50_long)
  df_D_post_long <- merge(df_D_post_long, df_D_post_q975_long)
  df_D_post_long <- merge(df_D_post_long, df_D_post_map_long)

  #names(df_D_post_long)
  df_D_post_long$time <- df_D_post_long$v_times
  df_D_post_long$v_times <- NULL

  #show_col(hue_pal()(10))
  # "#F8766D" "#7CAE00" "#00BFC4" "#C77CFF"
  v_col <- hue_pal()(10) # number of data sources plus MAP
  v_col[10] <- "#000000"

  colour_scale <- scale_colour_manual(name="", 
    values=c("Maximum a posterior" = v_col[10],
    "AgCensus"            = v_col[1], 
    "CORINE"              = v_col[2],
    "CROME"               = v_col[3],
    "CS"                  = v_col[4],
    "FC"                  = v_col[5],
    "IACS"                = v_col[6],
    "LCC"                 = v_col[7],
    "LCM"                 = v_col[8],
    "MODIS"               = v_col[9]))
  fill_scale <- scale_fill_manual(name="", 
    values=c("95% CI"              = v_col[10],
    "AgCensus"            = v_col[1], 
    "CORINE"              = v_col[2],
    "CROME"               = v_col[3],
    "CS"                  = v_col[4],
    "FC"                  = v_col[5],
    "IACS"                = v_col[6],
    "LCC"                 = v_col[7],
    "LCM"                 = v_col[8],
    "MODIS"               = v_col[9]))

  # plot D
  p <- ggplot(subset(df_D_post_long, time >= fig_start_time & time < fig_end_time), 
    aes(time, area_q50)) + theme_bw()
  p <- p + colour_scale
  p <- p + fill_scale
  # the interpolated observations
  p <- p + geom_ribbon(data = obs_unc$dt_D[time >= fig_start_time & time < fig_end_time],     
    aes(y    = area*(1-Fp), 
    ymin = area*(1-Fp) - sigma/5, 
    ymax = area*(1-Fp) + sigma/5, fill = data_source), alpha = 0.5)
  p <- p + geom_line  (data = obs_unc$dt_D[time >= fig_start_time & time < fig_end_time],     
    aes(y = area*(1-Fp), colour = data_source))
  # the UNinterpolated observations    
  p <- p + geom_point (data = obs_exc$dt_D[time >= fig_start_time & time < fig_end_time], 
    aes(y = area, colour = data_source))
  # the predictions
  p <- p + geom_ribbon(aes(ymin = area_q025, ymax = area_q975,  fill = "95% CI"), 
    alpha = 0.4)
  p <- p + geom_line(aes(time, area_MAP, colour = "Maximum a posterior"))
  p <- p + ylab(expression(paste(Area~km^2)))
  p <- p + ggtitle("Net Change")
  p <- p + facet_wrap(~ u, scales = "free_y")
  #p <- p + ylim(NA, 1000)
  p
  p_D <- p
  #ggsave(p, file = paste0(dir_output, "D_ts_UK.png"))

  #<!--- { DA2plotG -->
  # quantiles of area prediction from gross change
  df_G_post_q025 <- data.frame(v_times, t(a_G_post_q[1,,]))
  df_G_post_q50  <- data.frame(v_times, t(a_G_post_q[2,,]))
  df_G_post_q975 <- data.frame(v_times, t(a_G_post_q[3,,]))
  df_G_post_map  <- data.frame(v_times, t(a_G_post_map))
                        # melt(dt_G, id=c("time", "dtime"), variable.name = "u", value.name = "area")
  df_G_post_q025_long <- melt(df_G_post_q025, id=c("v_times"), variable.name = "u")
  df_G_post_q50_long  <- melt(df_G_post_q50 , id=c("v_times"), variable.name = "u")
  df_G_post_q975_long <- melt(df_G_post_q975, id=c("v_times"), variable.name = "u")
  df_G_post_map_long  <- melt(df_G_post_map,  id=c("v_times"), variable.name = "u")

  names(df_G_post_q025_long)[3] <- "area_q025"
  names(df_G_post_q50_long)[3]  <- "area_q50"
  names(df_G_post_q975_long)[3] <- "area_q975"
  names(df_G_post_map_long)[3]  <- "area_MAP"

  df_G_post_long <- merge(df_G_post_q025_long, df_G_post_q50_long)
  df_G_post_long <- merge(df_G_post_long, df_G_post_q975_long)
  df_G_post_long <- merge(df_G_post_long, df_G_post_map_long)

  df_G_post_long$time <- df_G_post_long$v_times
  df_G_post_long$v_times <- NULL

  #### new version
  # plot G
  p <- ggplot(subset(df_G_post_long, time >= fig_start_time & time < fig_end_time), aes(time, area_q50)) + theme_bw()
  p <- p + colour_scale
  p <- p + fill_scale
  # the interpolated observations    
  p <- p + geom_ribbon(data = obs_unc$dt_G[time >= fig_start_time & time < fig_end_time],     
    aes(y = area*(1-Fp),         
     ymin = area*(1-Fp) - sigma/5,         
     ymax = area*(1-Fp) + sigma/5, fill = data_source), alpha = 0.3)
  p <- p + geom_line(data = obs_unc$dt_G[time >= fig_start_time & time < fig_end_time],     
    aes(y = area*(1-Fp), colour = data_source))
  # the UNinterpolated observations    
  p <- p + geom_point (data = obs_exc$dt_G[time >= fig_start_time & time < fig_end_time], 
    aes(y = area, colour = data_source))
  # the predictions
  p <- p + geom_ribbon(aes(ymin = area_q025, ymax = area_q975,  fill = "95% CI"), 
    alpha = 0.4)
  p <- p + geom_line(aes(time, area_MAP, colour = "Maximum a posterior"))
  p <- p + ylab(expression(paste(Area~km^2)))
  p <- p + ggtitle("Gross Gain")
  p <- p + facet_wrap(~ u, scales = "free_y")
  p
  p_G <- p
  
  #<!--- { DA2plotL -->
  # quantiles of area prediction from gross change
  df_L_post_q025 <- data.frame(v_times, t(a_L_post_q[1,,]))
  df_L_post_q50  <- data.frame(v_times, t(a_L_post_q[2,,]))
  df_L_post_q975 <- data.frame(v_times, t(a_L_post_q[3,,]))
  df_L_post_map  <- data.frame(v_times, t(a_L_post_map))
                        # melt(dt_G, id=c("time", "dtime"), variable.name = "u", value.name = "area")
  df_L_post_q025_long <- melt(df_L_post_q025, id=c("v_times"), variable.name = "u")
  df_L_post_q50_long  <- melt(df_L_post_q50 , id=c("v_times"), variable.name = "u")
  df_L_post_q975_long <- melt(df_L_post_q975, id=c("v_times"), variable.name = "u")
  df_L_post_map_long  <- melt(df_L_post_map,  id=c("v_times"), variable.name = "u")

  names(df_L_post_q025_long)[3] <- "area_q025"
  names(df_L_post_q50_long)[3]  <- "area_q50"
  names(df_L_post_q975_long)[3] <- "area_q975"
  names(df_L_post_map_long)[3]  <- "area_MAP"

  df_L_post_long <- merge(df_L_post_q025_long, df_L_post_q50_long)
  df_L_post_long <- merge(df_L_post_long, df_L_post_q975_long)
  df_L_post_long <- merge(df_L_post_long, df_L_post_map_long)

  df_L_post_long$time <- df_L_post_long$v_times
  df_L_post_long$v_times <- NULL

  #### new version
  # plot L
  p <- ggplot(subset(df_L_post_long, time >= fig_start_time & time < fig_end_time), aes(time, area_q50)) + theme_bw()
  p <- p + colour_scale
  p <- p + fill_scale
  # the interpolated observations    
  p <- p + geom_ribbon(data = obs_unc$dt_L[time >= fig_start_time & time < fig_end_time],     
    aes(y  = area*(1-Fp),         
      ymin = area*(1-Fp) - sigma/5,         
      ymax = area*(1-Fp) + sigma/5, fill = data_source), alpha = 0.3)
  p <- p + geom_line(data = obs_unc$dt_L[time >= fig_start_time & time < fig_end_time],     
    aes(y = area*(1-Fp), colour = data_source))
  # the UNinterpolated observations    
  p <- p + geom_point (data = obs_exc$dt_L[time >= fig_start_time & time < fig_end_time], 
    aes(y = area, colour = data_source))
  # the predictions
  p <- p + geom_ribbon(aes(ymin = area_q025, ymax = area_q975,  fill = "95% CI"), 
    alpha = 0.4)
  p <- p + geom_line(aes(time, area_MAP, colour = "Maximum a posterior"))
  p <- p + ylab(expression(paste(Area~km^2)))
  p <- p + ggtitle("Gross Loss")
  p <- p + facet_wrap(~ u, scales = "free_y")
  p
  p_L <- p
  
  #<!--- { DA2plotB -->
  # quantiles of area prediction from gross change
  df_B_post_q025 <- data.frame(v_times, t(a_B_post_q[1,,]))
  df_B_post_q50  <- data.frame(v_times, t(a_B_post_q[2,,]))
  df_B_post_q975 <- data.frame(v_times, t(a_B_post_q[3,,]))
  df_B_post_map  <- data.frame(v_times, t(a_B_post_map))
  #df_B_post_map <- t(apply(a_B_post[,,], FUN = as.vector, MARGIN = 3))

                        # melt(dt_G, id=c("time", "dtime"), variable.name = "u", value.name = "area")
  df_B_post_q025_long <- melt(df_B_post_q025, id=c("v_times"), variable.name = "u")
  df_B_post_q50_long  <- melt(df_B_post_q50 , id=c("v_times"), variable.name = "u")
  df_B_post_q975_long <- melt(df_B_post_q975, id=c("v_times"), variable.name = "u")
  df_B_post_map_long  <- melt(df_B_post_map,  id=c("v_times"), variable.name = "u")

  names(df_B_post_q025_long)[3] <- "area_q025"
  names(df_B_post_q50_long)[3]  <- "area_q50"
  names(df_B_post_q975_long)[3] <- "area_q975"
  names(df_B_post_map_long)[3]  <- "area_MAP"

  df_B_post_long <- merge(df_B_post_q025_long, df_B_post_q50_long)
  df_B_post_long <- merge(df_B_post_long, df_B_post_q975_long)
  df_B_post_long <- merge(df_B_post_long, df_B_post_map_long)

  # get the unchanging areas on the diagonal and set to NA
  names_diag <- names_u_matrix[c(1,8,15,22,29,36)]
  ind <- df_B_post_long$u %in% names_diag
  df_B_post_long[ind,] <- NA
  df_B_post_long <- na.omit(df_B_post_long)

  # split the land use change in to from and to
  fromTo <- reshape2::colsplit(string= as.character(df_B_post_long$u), 
    pattern="_to_", names=c("u_from", "u_to"))
  df_B_post_long <- data.frame(df_B_post_long, fromTo)

  df_B_post_long <- within(df_B_post_long, 
    land_use_change <- factor(u, 
    levels = names_u_matrix_byrow))
  df_B_post_long <- within(df_B_post_long, 
    u_from <- factor(u_from, 
    levels = names_u))
  df_B_post_long <- within(df_B_post_long, 
    u_to <- factor(u_to, 
    levels = names_u))
  #with(df_B_post_long, levels(u_from))

  df_B_post_long$time <- df_B_post_long$v_times
  df_B_post_long$v_times <- NULL

  # p <- ggplot(df_B_post_long, aes(time, area_q50))
  # p <- p + colour_scale
  # p <- p + fill_scale
  # p <- p + geom_ribbon(aes(ymin = area_q025, ymax = area_q975, 
       # fill = "95% CI"), alpha = 0.3)
  # # p <- p + geom_point(data = obs_unc$dt_B[time >= 1990 & time < 2020], 
    # # aes(y = area*(1-Fp), colour = data_source))
  # p <- p + geom_line(aes(time, area_MAP, colour = "Maximum a posterior"))
  # p <- p + ylab(expression(paste(Area~km^2))) + xlab("Year")
  # p <- p + facet_grid(u_from ~ u_to, scales = "free_y")
  # p
  #ggsave(p, file = paste0(runDir, "beluc_B_UK.png"))

  #### new version
  # plot B
  p <- ggplot(subset(df_B_post_long, time >= fig_start_time & time < fig_end_time), aes(time, area_q50)) + theme_bw()
  p <- p + colour_scale
  p <- p + fill_scale
  # the interpolated observations    
  p <- p + geom_ribbon(data = obs_unc$dt_B[time >= fig_start_time & time < fig_end_time],     
    aes(y  = area*(1-Fp),         
      ymin = area*(1-Fp) - sigma/50,         
      ymax = area*(1-Fp) + sigma/50, fill = data_source), alpha = 0.3)
  p <- p + geom_line(data = obs_unc$dt_B[time >= fig_start_time & time < fig_end_time],     
    aes(y = area*(1-Fp), colour = data_source))
  # the UNinterpolated observations    
  p <- p + geom_point (data = obs_exc$dt_B[time >= fig_start_time & time < fig_end_time], 
    aes(y = area, colour = data_source))
  # the predictions
  p <- p + geom_ribbon(aes(ymin = area_q025, ymax = area_q975,  fill = "95% CI"), 
    alpha = 0.4)
  p <- p + geom_line(aes(time, area_MAP, colour = "Maximum a posterior"))
  p <- p + ylab(expression(paste(Area~km^2)))
  p <- p + ggtitle("Beta matrix")
  p <- p + facet_grid(u_from ~ u_to, scales = "free_y")
  p
  p_B <- p
  
  pdf(file = paste0(dir_output, "/post_plots_UK.pdf"))
   print(p_B)
   print(p_G)
   print(p_L)
   print(p_D)
  dev.off()

  saveRDS(  post$df_B, file = paste0(dir_output, "/df_B_post.rds"))
  saveRDS(  post$df_G, file = paste0(dir_output, "/df_G_post.rds"))
  saveRDS(  post$df_L, file = paste0(dir_output, "/df_L_post.rds"))
  saveRDS(  post$df_D, file = paste0(dir_output, "/df_D_post.rds"))
  write.csv(post$df_B, file = paste0(dir_output, "/df_B_post.csv"), row.names = FALSE)
  write.csv(post$df_G, file = paste0(dir_output, "/df_G_post.csv"), row.names = FALSE)
  write.csv(post$df_L, file = paste0(dir_output, "/df_L_post.csv"), row.names = FALSE)
  write.csv(post$df_D, file = paste0(dir_output, "/df_D_post.csv"), row.names = FALSE)

  return(list(p_B = p_B, p_G = p_G, p_L = p_L, p_D = p_D,
    df_B = df_B_post_long,
    df_G = df_G_post_long,
    df_L = df_L_post_long,
    df_D = df_D_post_long
    ))
}
