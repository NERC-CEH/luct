## ----luc_track_pkg, eval=TRUE------------------------------------------------
#' Processing and analysis of Land Use Vectors
#'
#' luc_track provides spatial utility functions and data.
"_PACKAGE"
#> [1] "_PACKAGE"

# this will load a host of packages used across processing scripts, in one location.

# pkgs <- c("units", "tidyr", "sp", "sf", "raster", "rgeos", "rgdal", "grid", "data.table", "spCEH", "scico", "ggplot2", "ggforce", "plyr", "dplyr", "stars", "BayesianTools", "ggthemes")
# pkgs_available <- sapply(pkgs, require, character.only = TRUE)
# if (!all(pkgs_available)) stop(paste("Not all packages available - check: ", which(!pkgs_available)))

#' Function to calculate the time series of areas of each land use
#' from a stars object or RasterStack converted to a data table
#' Get the area time series from a data table.
#' 
#' @param dt_U A land use stack.
#' @param cellArea Grid cell area, preferably with units attribute.
#' @return A matrix of area by time (rows) and land use (cols).
#' @export
#' @examples
#' dt_A <- getAreaTimeSeries(dt_U, cellArea)
getAreaTimeSeries <- function(dt_U, cellArea){
  dt_A <- dt_U[, .N, by = c("time", "u")]
  dt_A$area <- dt_A$N * cellArea
  return(dt_A)
}

#' Function to calculate the net change in area 
#'  for each land use from a land-use transition matrix
#'
#' @param v_B A land-use transition matrix, or its vector form
#' @param n_u The number of land-use classes
#' @return A vector of the net change in area by land use
#' @export
#' @examples
#' df_A <- getAreaTimeSeries_fromBetaArray(a_B, names_u, v_times)
getAreaTimeSeries_fromBetaArray <- function(a_B, 
  names_u = names_u, v_times = 1:dim(a_B)[3]){
  n_u <- dim(a_B)[1]
  n_t <- dim(a_B)[3]
  if (n_u != length(names_u)) stop("Mismatch in number of land-use classes.")
  if (n_t != length(v_times)) stop("Mismatch in number of time steps.")
  df_A <- data.frame(matrix(0, nrow = n_t, ncol = n_u))
  names(df_A) <- names_u
  # rowSums gives area at t-1
  # so for first time
  df_A[1,] <- rowSums(a_B[,, 2]) 
  for (i_t in 2:n_t){
    df_A[i_t,] <- colSums(a_B[,, i_t])
  }
  df_A$time <- v_times
  # do we want a data frame or a data table, or just a matrix?
  dt_A <- data.table(df_A)
  return(dt_A)  
}

#' Function to calculate the transition matrices
#' from a stars object or RasterStack converted to a data table
#' Get the Beta matrices series from a data table.
#' 
#' @param dt_U A data table containing the land use at every location.
#' @param cellArea Area of single cell in the raster grid, preferably with units attribute, which defines output area unit.
#' @param removeDiag Remove data for unchanging land use, so set the matrix diagonal to zeroes.
#' @return An array of matrices of area by time (rows) and land use (cols). If rateOfChange is TRUE, units are area/year, despite the attribute saying area
#' @export
#' @examples
#' a_B <- getBetaMatrices(dt_U, cellArea)
getBetaMatrices <- function(dt_U, names_u = names_u, cellArea, 
  removeDiag = FALSE, rateOfChange = TRUE){
  v_times <- unique(dt_U$time)
  n_u <- length(names_u)
  n_t <- length(v_times)
  # we need the lu classes to be integers for this algorithm
  if (!is.numeric(dt_U$u)) dt_U$u  <- match(dt_U$u, names_u) 
  v_dtime <- c(0, diff(v_times))
  a_B_ijt <- array(0, c(n_u, n_u, n_t))

  for (i_t in 2:n_t){
    i_time     <- v_times[i_t]
    i_timeprev <- v_times[(i_t-1)]
    v_U_t     <- dt_U[time == i_time, u]
    v_U_tprev <- dt_U[time == i_timeprev, u]
    for (i_u in 1:n_u) {
      for (j_u in 1:n_u) {
        a_B_ijt[i_u,j_u,i_t] <- 
          sum(v_U_tprev == i_u 
            & v_U_t     == j_u, na.rm = TRUE)
      } # j_u
    }   # i_u
    if (removeDiag) diag(a_B_ijt[,,i_t]) <- 0
    if (rateOfChange){
      # divide by time increment, assumed years
      m_dtime <- v_dtime[i_t] + diag(n_u)
      # but not diagonal - unchanging elements
      diag(m_dtime) <- 1
      a_B_ijt[,,i_t] <- a_B_ijt[,,i_t] / m_dtime
    } # rateOfChange
  }   # itime
  a_B_ijt <- a_B_ijt * cellArea # units of cellArea per year
  return(a_B_ijt)
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

#' Function to calculate the time series of change for each land use
#' from an array of Beta matrices. Change is returned in the form of
#' Beta matrices, gross Losses, net change in Area, and gross Gains,
#' collectively "BLAG".
#' 
#' @param a_B An array of Beta matrices containing areas changing between intervals t and t-1.
#' @param names_u A vector of the names of land-use classes.
#' @param v_times A vector of the times corresponding to the time dimension in a_B, assumed to be years.
#' @return A BLAG object.
#' @export
#' @examples
#' blag <- getBLAG(a_B, names_u)
getBLAG <- function(a_B, names_u = names_u, name_data_source = "",
      region = NA, v_times = 1:dim(a_B)[3], cellArea = cellArea){
  n_u <- dim(a_B)[1]
  n_t <- dim(a_B)[3]
  if (n_u != length(names_u)) stop("Mismatch in number of land-use classes.")
  if (n_t != length(v_times)) stop("Mismatch in number of time steps.")

  # Need B matrix in rate units  (e.g. m^2/y) for comparisons where time intervals are not annual
  # Elsewhere need in area units (e.g. m^2)
  a_B_rate <- getBetaRateMatrices(a_B, names_u = names_u, 
    v_times = v_times, removeDiag = TRUE)
  #** this is not giving rates
  # B data table
  dt_B <- as.data.table(remove_units(a_B_rate)) # first time - no difference to calculate, so all zeroes
  names(dt_B) <- c("u_from", "u_to", "time", "area")
  dt_B$area <- remove_units(dt_B$area)

  dt_B$time <- v_times[dt_B$time]
  dt_B$u_from <- as.factor(dt_B$u_from)
  dt_B$u_to   <- as.factor(dt_B$u_to)
  levels(dt_B$u_from) = names_u
  levels(dt_B$u_to) = names_u

  # A time series
  dt_A_wide <- getAreaTimeSeries_fromBetaArray(a_B, names_u, v_times)
  dt_A <- dt_A_wide %>% gather(u, area, one_of(names_u))
  dt_A <- na.omit(dt_A, cols="u")
  dt_A$u <- factor(dt_A$u, levels = names_u)

  # D time series
  dt_D <- as.data.table(t(apply(a_B[,,], MARGIN = 3, FUN = getAreaNetChange_fromBeta, n_u = n_u)))
  names(dt_D) <- names_u
  dt_D <- data.table(time = v_times, dtime = c(0, diff(v_times)), dt_D)
  dt_D <- melt(dt_D, id=c("time", "dtime"), variable.name = "u", value.name = "area")
  dt_D$area <- dt_D$area / dt_D$dtime # convert km2 to km2/year
  dt_D$dtime <- NULL
  
  # G time series
  # calc gross gain from B; MARGIN = 3 because that is the time dimension
  dt_G <- as.data.table(t(apply(a_B[,,], MARGIN = 3, FUN = getAreaGrossGain_fromBeta, n_u = n_u)))
  names(dt_G) <- names_u
  dt_G <- data.table(time = v_times, dtime = c(0, diff(v_times)), dt_G)
  dt_G <- melt(dt_G, id=c("time", "dtime"), variable.name = "u", value.name = "area")
  dt_G$area <- dt_G$area / dt_G$dtime # convert km2 to km2/year
  dt_G$dtime <- NULL
  
  # L time series
  # calc gross loss from B; MARGIN = 3 because that is the time dimension
  dt_L <- as.data.table(t(apply(a_B[,,], MARGIN = 3, FUN = getAreaGrossLoss_fromBeta, n_u= n_u)))
  names(dt_L) <- names_u
  dtime = c(0, diff(v_times))
  dt_L <- data.table(time = v_times, dtime = c(0, diff(v_times)), dt_L)
  dt_L <- melt(dt_L, id=c("time", "dtime"), variable.name = "u", value.name = "area")
  dt_L$area <- dt_L$area / dt_L$dtime # convert km2 to km2/year
  dt_L$dtime <- NULL
  
  # # D, net change in area
  # # calculate annual change by difference
  # dt_D_wide <- as.data.table(apply(dt_A_wide, MARGIN = 2, diff))
  # # diff produces n-1 rows, so need to add a row of zeroes to conform
  # dt_D_wide <- rbind(dt_D_wide[1,], dt_D_wide)
  # dt_D_wide[1,] <- 0
  # dt_D_wide$dtime <- dt_D_wide$time
  # dt_D_wide$time  <- v_times
  
  # dt_D <- melt(dt_D_wide, id=c("time", "dtime"), variable.name = "u", value.name = "area")
  # dt_D$area <- dt_D$area / dt_D$dtime # convert km2 to km2/year
  # dt_D$dtime <- NULL

  # reshape to wide format
  # G and L, Gross Gains and Losses
  dt_G_wide <- dcast(dt_G, time ~ u, value.var="area")
  rownames(dt_G_wide) <- dt_G_wide$time
  dt_G_wide$time <- NULL
  dt_L_wide <- dcast(dt_L, time ~ u, value.var="area")
  rownames(dt_L_wide) <- dt_L_wide$time
  dt_L_wide$time <- NULL
  # D
  dt_D_wide <- dcast(dt_D, time ~ u, value.var="area")
  rownames(dt_D_wide) <- dt_D_wide$time
  dt_D_wide$time <- NULL

  # add data source name as a column
  dt_A$data_source <- name_data_source
  dt_D$data_source <- name_data_source
  dt_B$data_source <- name_data_source
  dt_G$data_source <- name_data_source
  dt_L$data_source <- name_data_source
  
  # add region name as a column
  dt_A$region <- region
  dt_D$region <- region
  dt_B$region <- region
  dt_G$region <- region
  dt_L$region <- region

  return(list(
    v_times = v_times,
    v_dtime = dtime,
    dt_A = dt_A, 
    dt_D = dt_D, 
    dt_B = dt_B, 
    a_B = a_B, 
    dt_G = dt_G, 
    dt_L = dt_L,
    # reshaped
    dt_A_wide  = dt_A_wide ,
    dt_D_wide = dt_D_wide,
    dt_G_wide  = dt_G_wide ,
    dt_L_wide  = dt_L_wide
  ))
}

#' Function to calculate the Beta matrices of land-use change
#' in the form of rates, i.e. area per year, from an array of Beta matrices
#' and a time vector. This is just divides by the time difference, 
#' so only matters when intervals are irregular or not annual. 
#' This can surely be done more neatly with matrix division.
#' 
#' @param a_B An array of Beta matrices containing areas changing between intervals t and t-1.
#' @param names_u A vector of the names of land-use classes.
#' @param v_times A vector of the times corresponding to the time dimension in a_B, assumed to be years.
#' @param removeDiag Logical - set matrix diagonal (unchanging area to zero).
#' @return An array of Beta matrices containing the rate of change (area/year) changing between intervals t and t-1.
#' @export
#' @examples
#' a_B_rate <- getBetaRateMatrices(a_B, names_u)
getBetaRateMatrices <- function(a_B, names_u = names_u,
  v_times = 1:dim(a_B)[3], removeDiag = FALSE){
  n_u <- dim(a_B)[1]
  n_t <- dim(a_B)[3]
  if (n_u != length(names_u)) stop("Mismatch in number of land-use classes.")
  if (n_t != length(v_times)) stop("Mismatch in number of time steps.")
  v_dtime <- c(0, diff(v_times))
  a_B_ijt <- a_B

  for (i_t in 2:n_t){
    if (removeDiag) diag(a_B_ijt[,,i_t]) <- 0
    # divide by time increment, assumed years
    m_dtime <- v_dtime[i_t] + diag(n_u)
    # but not diagonal - unchanging elements
    diag(m_dtime) <- 1
    a_B_ijt[,,i_t] <- a_B_ijt[,,i_t] / m_dtime
  }   # itime
  return(a_B_ijt)
}

#' Function to create data structure with all terms 
#' Beta matrix, Losses, Area time series, and Gains
#' (B, L, A, G) from U. U is a stars object with integer 
#' land-use classes (u) in space and time dimensions.

#' @param   v_times Numeric vector of times (in years).
#' @param   v_fnames Character vector of file names.
#' @param   name_data_source Character string for name of data source.
#' @param   names_u Character vector of names of land uses.
#' @param   returnU Logical Whether to return the whole U object instead. Default is FALSE.
#' @return A BLAG object.
#' @export
#' @examples
#' blag <- getBLAG_fromU(v_times = c(2000, 2006),
#'   v_fnames = c("data/CORINE/Level1/r_U_cor_1000m_2000.tif", 
#'        "data/CORINE/Level1/r_U_cor_1000m_2006.tif"),
#'   name_data_source = "CORINE", 
#'   region = "sc",
#'   names_u = names_u)
getBLAG_fromU <- function(
  v_times,
  v_fnames,
  name_data_source, 
  region = "uk",
  names_u, 
  returnU = FALSE){
  
  file.exists(v_fnames)
  s_U <- stack(v_fnames)
  st_U <- st_as_stars(s_U)
  st_U <- st_set_dimensions(st_U, "band", values = v_times, names = "time")
 
  n_u <- length(names_u)
  n_t <- dim(st_U)[3]
  v_times <- st_get_dimension_values(st_U, 3)
  if (names(dim(st_U)[3]) != "time") stop("Third dimension of stars object should be 'time'.")  
  # in which case need to correct this using line below:
  #st_U <- st_set_dimensions(st_U, "band", values = v_times, names = "time")
  cellArea <- mean(st_area(st_U)$area)
  dt_U <- as.data.table(st_U)
  names(dt_U)[4] <- "u"
  # convert indices to names, in case any are missing
  dt_U$u <- names_u[dt_U$u] 
  levels(dt_U$u) <- names_u

  # mask out regions to be excluded
  # read the mask file with the matching resolution
  # region is called "country" in mask file - to be fixed
  fname <- paste0("./data-raw/mask/dt_land_", res(s_U)[1], "m.qs")
  dt_mask <- qread(fname)
  v_region_id <-  c(1, 2, 3, 4)
  v_region_all <-  c("en", "sc", "wa", "ni", "uk")
  country_id <- match(region, v_region_all)
  if (region == "uk") country_id <- v_region_id
  dt <- merge(dt_mask, dt_U)
  dt <- dt[country %in% country_id]
  # dim(dt)
  # names(dt) # "x"    "y"    "time" "u"
  dt_U <- dt[, .(x, y, time, u)]
  
  # B matrix
  a_B <- getBetaMatrices(dt_U, names_u = names_u, cellArea = cellArea, 
    removeDiag = FALSE, rateOfChange = FALSE)
  # a_B_rate <- getBetaMatrices(dt_U, names_u = names_u, cellArea = cellArea, 
    # removeDiag = TRUE, rateOfChange = TRUE)
  # a_B <- remove_units(a_B)
  
  #** get this correct
  blag <- getBLAG(a_B, names_u = names_u, name_data_source = name_data_source,
    region = region, v_times = v_times, cellArea = cellArea)

  # make list of items to return
  l <- list(
    cellArea = cellArea, 
    v_times = blag$v_times,
    v_dtime = blag$v_dtime,
    dt_A  = blag$dt_A, 
    dt_D = blag$dt_D,
    dt_B  = blag$dt_B, 
    a_B   = a_B, 
    dt_G  = blag$dt_G, 
    dt_L  = blag$dt_L,
    # reshaped
    dt_A_wide  = blag$dt_A_wide ,
    dt_D_wide = blag$dt_D_wide,
    dt_G_wide  = blag$dt_G_wide ,
    dt_L_wide  = blag$dt_L_wide
  )
  # do not usually want to return this - too big
  if (returnU) {
    l <- list(l, 
      st_U = st_U, 
      dt_U = dt_U)
  }
  return(l)
}

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


#' Function to create a diverging colour palette for a RasterLayer
#'   together with raster attribute table and time series data frame
#'
#' @param r A RasterLayer
#' @param start.color Starting colour
#' @param end.color   Ending colour
#' @param mid.value   Mid-range value
#' @param mid.color   Mid-range colour
#' @return A colour palette
#' @export
#' @examples
#' col <- diverge.color(r)
diverge.color <- function(r, start.color = "red", end.color = "blue", mid.value=0, mid.color="light grey"){
  # based on ideas from Maureen Kennedy, Nick Povak, and Alina Cansler

  # creates a palette for the current session for a divergent-color
  # graphic with a non-symmetric range
  # "cuts" = the number of slices to be made in the range above and below "mid.value"
  # vector version is quicker than raster::cellStats
  v <- getValues(r)
  min.value <- min(v, na.rm=TRUE)
  max.value <- max(v, na.rm=TRUE)
  #min.value <- cellStats(r, "min")
  #max.value <- cellStats(r, "max")

  ramp1 <- colorRampPalette(c(start.color,mid.color))
  ramp2 <- colorRampPalette(c(mid.color,end.color))

  # now specify the number of values on either side of "mid.value"

  max.breaks <- round(max.value - mid.value)
  max.breaks <- max(max.breaks, 5) # in case rounding removes too many breaks
  min.breaks <- round(mid.value - min.value)
  min.breaks <- max(min.breaks, 5) # in case rounding removes too many breaks

  num.breaks <- max(max.breaks,min.breaks)

  low.ramp <- ramp1(num.breaks)
  high.ramp <- ramp2(num.breaks)

  # now create a combined ramp from the higher values of "low.ramp" and 
  # the lower values of "high.ramp", with the longer one using all values 
  # high.ramp starts at 2 to avoid duplicating zero

  myColors <- c(low.ramp[(num.breaks-min.breaks):num.breaks],high.ramp[2:max.breaks])
  return(myColors)
}

#' Function to derive a list of all unique vectors and their areas from the b_U brick
#' Get the unique vectors from a b_U brick.
#' 
#' @param b_U A RasterBrick of land use U.
#' @return A vector object containing a data frame and a raster.
#' @export
#' @examples
#' v <- getVectors(b_U)
#' v <- clusterR(b_U, getVectors)
getVectors <- function(b_U){
# b_U <- lu_iacs
  # change "NA" to "0" - needs to be single char for this to work
  naToZero <- function(x) { x[is.na(x)] <- 0; return(x)} 
  #freq(b_U)
  b_U <- calc(b_U, naToZero)
  
  # concat each year's b_U integer code as a char string eg. 111222333
  # named luv_int: land use vector represented as integers
  ## rename luv_int_ch: char representation of integer - it exceeds size of max.integer
  ## declare size of luv_int to speed it up?
  luv_int <- NULL
  for (iyr in 1:nlayers(b_U)){
    v <- getValues(b_U[[iyr]])
    luv_int <- paste(luv_int, v, sep="")
  }


  # find vectors with any missing values (coded as zero) and set to single zero
  containsZero <- grepl("0", luv_int)
  luv_int[containsZero] <- "0"

  r_luv <- setValues(b_U[[1]], as.factor(luv_int))

  ## make a data frame with the vectors put into cols
  lu_vectors <- table(luv_int)
  lu_vectors_df <- as.data.frame(lu_vectors)
  # find and remove any vectors with missing values (coded as zero)
  containsZero <- grepl("0", lu_vectors_df[,1])
  lu_vectors_df <- lu_vectors_df[!containsZero,]
  names(lu_vectors_df)
  lu_fac <- lu_vectors_df[,1]
  lu_vectors_split <- with(lu_vectors_df, strsplit(as.character(lu_vectors_df$luv_int), split = character(0))) 
  lu_vectors_split <- as.data.frame(lu_vectors_split)
  names(lu_vectors_split) <- NULL # names too long for file i/o
  lu_vectors_split <- t(lu_vectors_split)
  # vector area in km2 is number of 100x100 m cells / 100 -> km2
  ## this line previously assumed fixed resolution of 100 m 
  # should make area_gridcell_m2 = f(resolution)
  area_lu_vector <- lu_vectors_df$Freq * area_gridcell_m2 * m2ToKm2
  lu_df <- cbind.data.frame(lu_fac, area_lu_vector, lu_vectors_split)
  return(list(lu_df=lu_df, r_luv=r_luv))
}

names_u_chess <- c("white", "empty", "black")
names_u  <- c("woods", "crops", "grass", "rough", "urban", "other")
n_u  <- length(names_u)
colour_u <- c("lightblue", "green", "pink", "purple", "orange", "midnightblue")
