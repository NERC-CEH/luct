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

## ---- wrangle_AgCensus_Eng

#' Function to wrangle AgCensus data 
#'  from text files to R objects
#'
#' @param v_fpath Vector of filepaths to data files
#' @return A BLAG object
#' @export
#' @examples
#' fpath1 = "./data-raw/AgCensus/England/AgCensus_England_ha_1900-2010.csv"
#' fpath2 = "./data-raw/AgCensus/England/AgCensus_England_ha_1983-2019.csv"
#' x <- wrangle_AgCensus_Eng(c(fpath1, fpath2))$df_A
wrangle_AgCensus_Eng <- function(v_fpath = 
  c("./data-raw/AgCensus/England/AgCensus_England_ha_1900-2010.csv",
    "./data-raw/AgCensus/England/AgCensus_England_ha_1983-2019.csv")){

#<!--- { read_data_England1 -->
  df <- read.csv(v_fpath[1], na.strings = c("n/c", "-"))

  df$woods <-NA
  df$crops <- # crops are in columns 3-28
    rowSums(df[, 3:28], na.rm = TRUE)
  df$grass <- 
    df$Permanent.grass..excl.rough.grazing.
  df$rough <- 
    df$Total.rough.grazing
  df <- df[, c("year", "woods", "crops", "grass", "rough" )]

  # interpolate between decadal points
  df_allyears <- data.frame(year = 1900:1990)
  df <- merge(df_allyears, df, all.x = TRUE)
  df[,-1] <- na.approx(df[,-1])

  # set units to hectares (displayed in Excel as kha, but stored as ha)
  df[,-1] <- df[,-1] %>% map2_dfc("ha", ~set_units(.x, .y, mode = "standard"))
  # and convert units to km2
  df[,-1] <- df[,-1] %>% map2_dfc("km^2", ~set_units(.x, .y, mode = "standard"))

  df_hist <- df # make a copy of original data

#<!--- { read_data_England2 -->
  df <- read.csv(v_fpath[2], na.strings = c("n/c"))
  #df[is.na(df)] <- 0 # otherwise totals become NA
  #df[df == 0] <- NA # otherwise totals become NA

  df$woods <- 
    df$Woodland
  df$crops <- 
    df$Total.crops + 
    df$Uncropped.arable.land..No.set.aside.1983.1989...f.
  # need more complex form for addition with NAs  
  df$grass <- rowSums(df[,c(
    "Temporary.grass..sown.in.the.last.5.years.", 
    "Land.used.for.outdoor.pigs..g.", 
    "Grass.over.5.years.old")], na.rm=TRUE)
  df$rough <- 
    df$Common.rough.grazing..e. + 
    df$Sole.right.rough.grazing
  df <- df[, c("year", "woods", "crops", "grass", "rough" )]

  # code the first year as a known change point (better than in the data file), but still ugly
  df$year[1] <- df$year[1] + 0.5

  # set units to hectares (displayed in Excel as kha, but stored as ha)
  df[,-1] <- df[,-1] %>% map2_dfc("ha", ~set_units(.x, .y, mode = "standard"))
  # and convert units to km2
  df[,-1] <- df[,-1] %>% map2_dfc("km^2", ~set_units(.x, .y, mode = "standard"))

  # remove overlap between hist and contemporary data, which starts in 1983
  # remove overlap with contemporary data
  df_hist <- subset(df_hist, year <= 1983)
  df <- rbind(df_hist, df)

  df_orig <- df # make a copy of original data
 
  #<!--- { remove_steps -->
  knownChangePoint <- c(1983, 2009)
  ## 1977 and 2006 appear to be step changes but there is no record of a change in method
  unknownChangePoint <- c(1999, 2000, 2001, 2006) # NULL or add apparent step changes here
  # need to repeat for all variables
  df <- drop_units(df) # gam doesn't work with units
  for (j in 2:(length(df))){   # for each variable except year
  #j = 2  
    # go through all the columns except first (= year)
    # first, adjust for step changes where values are given under old and new basis
    for (i in length(knownChangePoint):1){
    #i = 1
      df_pre  <- subset(df, year <= knownChangePoint[i])
      df_post <- subset(df, year > knownChangePoint[i])
      obs_pre  <- df[,j][df$year == knownChangePoint[i]]
      obs_post <- df[,j][df$year == (knownChangePoint[i]+0.5)]
      jump <- obs_post - obs_pre
      df_pre[,j] <- df_pre[,j] + jump
      df <- rbind(df_pre, df_post)
    } # for each knownChangePoint

    # gam version for predicting area in changepoint year, based on a gam fitted to the previous data
    # use this where there are step changes, but values not given on old and new basis
    if (length(unknownChangePoint) > 0){
      for (i in length(unknownChangePoint):1){
      #i = 2
        df_pre  <- subset(df, year < unknownChangePoint[i])
        df_post <- subset(df, year >= unknownChangePoint[i])
        if ( sum(!is.na(df_pre[,j])) > 10 ){ # only do if >10 values, excluding NAs
          #k <- max(length(df_pre[,1])-1, 4)
          k <- min(length(df_pre[,1])-1, 10)
          mod_pre  <- gam(df_pre[,j] ~ s(year, bs = "cr", k = k), data = df_pre) 
          #k <- max(length(df_post[,1])-1, 4)
          #k <- min(length(df_post[,1])-1, 10)
          #mod_post <- gam(grassArea ~ s(year, bs = "cr", k = k), data = df_post)
          pred_cpt_pre  <- predict.gam(mod_pre, data.frame(year=unknownChangePoint[i]))
          #pred_cpt_post <- predict.gam(mod_post, data.frame(year=unknownChangePoint[i]))
          #pred_pre  <- predict.gam(mod_pre, df)
          obs_pre  <- df[,j][df$year == unknownChangePoint[i]]
          jump <- obs_pre - pred_cpt_pre
          df_pre[,j] <- df_pre[,j] + jump
          df <- rbind(df_pre, df_post)
        } # >10 values
      }   # i change points
    }     # if there are any change points
  }       # for each variable except year

  # remove duplicated rows with knownChangePoints, where one has year = XXXX.5
  #df <- subset(df, year %% floor(year) == 0) # so modulo division is not zero

  # restore units
  df[,-1] <- df[,-1] %>% map2_dfc("km^2", ~set_units(.x, .y, mode = "standard"))
  df_long      <- melt(df, id=c("year"), variable.name = "land_use", value.name = "area")
  df_long_orig <- melt(df_orig, id=c("year"), variable.name = "land_use", value.name = "area")
  names(df_long_orig) <- c("time", "u", "area" )
  names(df_long) <- c("time", "u", "area" )
  notDuplictates <- (df_long$time - as.integer(df_long$time)) == 0
  df_long <- df_long[notDuplictates, ]
  df_long$country <- "E"
  df_A_Eng <- df_long
  return(list(df_A = df_A_Eng, df_A_orig = df_long_orig))
}

## ---- wrangle_AgCensus_Sco

#' Function to wrangle AgCensus data 
#'  from text files to R objects
#'
#' @param v_fpath Vector of filepaths to data files
#' @return A BLAG object
#' @export
#' @examples
#' fpath1 = "./data-raw/AgCensus/Scotland/AgCensus_Scotland_ha_1883-2014.csv"
#' fpath2 = "./data-raw/AgCensus/Scotland/AgCensus_Scotland_ha_2009-2019.csv"
#' x <- wrangle_AgCensus_Sco(c(fpath1, fpath2))
wrangle_AgCensus_Sco <- function(v_fpath = 
  c("./data-raw/AgCensus/Scotland/AgCensus_Scotland_ha_1883-2014.csv",
    "./data-raw/AgCensus/Scotland/AgCensus_Scotland_ha_2009-2019.csv")){

  # read historical data
  df <- read.csv(v_fpath[1], na.strings = c("n/a", ":", "b"))
  df_hist <- df # make a copy of original data

  # read current data
  df <- read.csv(v_fpath[2], na.strings = c("n/a", ":", "b"))

  df$woods <- 
    df$Other.land..including.woodland..3.
  df$crops <- 
    df$Total.crops..fallow..and.set.aside.
  df$grass <- 
    df$Total.grass
  df$rough <- 
    df$Rough.grazing
  df <- df[, c("year", "woods", "crops", "grass", "rough" )]

  # after 2014, woods area is not separated from other non-agric land included in the census.
  # as a crude fix, assume the non-woods area stays constant after 2014
  # df_hist$woods[df_hist$year == 2014] # woods area in 2014
  # df$woods[df$year == 2014]           # woods + other non-woods area in 2014
  area_NonWoods <- df$woods[df$year == 2014] - # the difference is the 
         df_hist$woods[df_hist$year == 2014]   # non-woods area
  df$woods <- df$woods - area_NonWoods  # subtract this from total woods + other non-woods reported after 2014

  # remove overlap between hist and contemporary data, which starts in 2009
  # the hist data are needed for the above comparison of wood/non wood area until this point
  df_hist <- subset(df_hist, year < 2009)
  df <- rbind(df_hist, df)

  # set units to hectares (displayed in Excel as kha, but stored as ha)
  df[,-1] <- df[,-1] %>% map2_dfc("ha",   ~set_units(.x, .y, mode = "standard"))
  # and convert units to km2
  df[,-1] <- df[,-1] %>% map2_dfc("km^2", ~set_units(.x, .y, mode = "standard"))

  df_orig <- df # make a copy of original data
  #df <- df_orig # make a copy of original data

  #<!--- { remove_steps_Scotland -->
  knownChangePoint <- c(1970, 1973, 1982)
  ## 1977 and 2006 appear to be step changes but there is no record of a change in method
  unknownChangePoint <- c(2009) # 2006, 1959, 1960, 1961, 1977, # NULL or add apparent step changes here
  # need to repeat for all variables
  df <- drop_units(df) # gam doesn't work with units
  for (j in 2:(length(df))){   # for each variable except year
  #j = 2  
    # go through all the columns except first (= year)
    # first, adjust for step changes where values are given under old and new basis
    for (i in length(knownChangePoint):1){
    #i = 1
      df_pre  <- subset(df, year <= knownChangePoint[i])
      df_post <- subset(df, year > knownChangePoint[i])
      obs_pre  <- df[,j][df$year == knownChangePoint[i]]
      obs_post <- df[,j][df$year == (knownChangePoint[i]+0.5)]
      jump <- obs_post - obs_pre
      df_pre[,j] <- df_pre[,j] + jump
      df <- rbind(df_pre, df_post)
    } # for each knownChangePoint

    # gam version for predicting area in changepoint year, based on a gam fitted to the previous data
    # use this where there are step changes, but values not given on old and new basis
    if (length(unknownChangePoint) > 0){
      for (i in length(unknownChangePoint):1){
      #i = 1
        df_pre  <- subset(df, year < unknownChangePoint[i])
        df_post <- subset(df, year >= unknownChangePoint[i])
        if ( sum(!is.na(df_pre[,j])) > 10 ){ # only do if >10 values, excluding NAs
          #k <- max(length(df_pre[,1])-1, 4)
          k <- min(length(df_pre[,1])-1, 10)
          mod_pre  <- gam(df_pre[,j] ~ s(year, bs = "cr", k = k), data = df_pre) 
          #k <- max(length(df_post[,1])-1, 4)
          #k <- min(length(df_post[,1])-1, 10)
          #mod_post <- gam(grassArea ~ s(year, bs = "cr", k = k), data = df_post)
          pred_cpt_pre  <- predict.gam(mod_pre, data.frame(year=unknownChangePoint[i]))
          #pred_cpt_post <- predict.gam(mod_post, data.frame(year=unknownChangePoint[i]))
          #pred_pre  <- predict.gam(mod_pre, df)
          obs_pre  <- df[,j][df$year == unknownChangePoint[i]]
          jump <- obs_pre - pred_cpt_pre
          df_pre[,j] <- df_pre[,j] + jump
          df <- rbind(df_pre, df_post)
        } # >10 values
      }   # i change points
    }     # if there are any change points
  }       # for each variable except year
  # remove duplicated rows with knownChangePoints, where one has year = XXXX.5
  #df <- subset(df, year %% floor(year) == 0) # so modulo division is not zero

  # restore units
  df[,-1] <- df[,-1] %>% map2_dfc("km^2", ~set_units(.x, .y, mode = "standard"))

  df_long      <- melt(df, id=c("year"), variable.name = "u", value.name = "area")
  df_long_orig <- melt(df_orig, id=c("year"), variable.name = "u", value.name = "area")

  names(df_long) <- c("time", "u", "area" )
  names(df_long_orig) <- c("time", "u", "area" )
  notDuplictates <- (df_long$time - as.integer(df_long$time)) == 0
  df_long <- df_long[notDuplictates, ]
  df_long$country <- "S"
  df_A_Sco <- df_long
  return(list(df_A = df_A_Sco, df_A_orig = df_long_orig))
}

## ---- wrangle_AgCensus_Wal

#' Function to wrangle AgCensus data 
#'  from text files to R objects
#'
#' @param v_fpath Vector of filepaths to data files
#' @return A BLAG object
#' @export
#' @examples
#' fpath1 = "./data-raw/AgCensus/Wales/AgCensus_Wales_ha_1867-2012.csv"
#' fpath2 = "./data-raw/AgCensus/Wales/AgCensus_Wales_ha_1998-2019.csv"
#' x <- wrangle_AgCensus_Wal(c(fpath1, fpath2))
wrangle_AgCensus_Wal <- function(v_fpath = 
  c("./data-raw/AgCensus/Wales/AgCensus_Wales_ha_1867-2012.csv",
    "./data-raw/AgCensus/Wales/AgCensus_Wales_ha_1998-2019.csv")){

  # read historical data
  df <- read.csv(v_fpath[1], na.strings = c("n/a", ":", "b"))

  df$woods <- NA
  df$crops <- 
    df$Crops
  df$grass <- 
    df$Permanent.pasture + 
    df$New.grass
  df$rough <- 
    df$Rough.grazing
  df <- df[, c("year", "woods", "crops", "grass", "rough" )]
  df_hist <- df # make a copy of original data

  # read current data
  df <- read.csv(v_fpath[2], na.strings = c("n/a"))
  # names(df)
  # head(df)
  # summary(df)
  #df[is.na(df)] <- 0 # otherwise totals become NA

  df$woods <- 
    df$Woodland
  df$crops <- rowSums(df[,c(
    "Section.sub.total", 
    "Section.sub.total.1", 
    "Set.aside...total")], na.rm=TRUE)
  df$grass <- 
    df$Grassland.under.5.years.old + 
    df$Grassland.over.5.years.old 
  df$rough <- 
    df$Common.rough.grazing + 
    df$Sole.rights.rough.grazing
  df <- df[, c("year", "woods", "crops", "grass", "rough" )]

  # remove overlap between hist and contemporary data, which starts in 1998
  df_hist <- subset(df_hist, year < 1998)
  df <- rbind(df_hist, df)

  # set units to hectares (displayed in Excel as kha, but stored as ha)
  df[,-1] <- df[,-1] %>% map2_dfc("ha",   ~set_units(.x, .y, mode = "standard"))
  # and convert units to km2
  df[,-1] <- df[,-1] %>% map2_dfc("km^2", ~set_units(.x, .y, mode = "standard"))

  df_orig <- df # make a copy of original data, before removing step changes

  df_long      <- melt(df, id=c("year"), variable.name = "land_use", value.name = "area")
  #df_long_orig <- melt(df_orig, id=c("year"), variable.name = "land_use", value.name = "area")

  names(df_long) <- c("time", "u", "area" )
  notDuplictates <- (df_long$time - as.integer(df_long$time)) == 0
  df_long <- df_long[notDuplictates, ]
  df_long$country <- "W"
  df_A_Wal <- df_long
  return(list(df_A = df_A_Wal, df_A_orig = df_A_Wal))
}

## ---- wrangle_AgCensus_NIr

#' Function to wrangle AgCensus data 
#'  from text files to R objects
#'
#' @param fpath Filepath to text file
#' @return A BLAG object
#' @export
#' @examples
#' fpath = "./data-raw/AgCensus/NIreland/AgCensus_NIreland_ha_1981-2019.csv"
#' x <- wrangle_AgCensus_NIr(fpath)
wrangle_AgCensus_NIr <- function(
  fpath = "./data-raw/AgCensus/NIreland/AgCensus_NIreland_ha_1981-2019.csv"){

  # read NIr data
  df <- read.csv(fpath, na.strings = c("n/a", ":", "b"))

  df$woods <- df$Woods.and.plantations
  df$crops <- 
    df$TOTAL.CROPS +
    df$Set.aside
  df$grass <- df$TOTAL.GRASS
  df$rough <- df$Hill.or.Rough.Land
  df <- df[, c("year", "woods", "crops", "grass", "rough" )]

  # set units to hectares
  df[,-1] <- df[,-1] %>% map2_dfc("ha",   ~set_units(.x, .y, mode = "standard"))
  # and convert units to km2
  df[,-1] <- df[,-1] %>% map2_dfc("km^2", ~set_units(.x, .y, mode = "standard"))
  #dim(df)
  #sum(df[39,-1]) # total 10129.21 km^2 cf. total area NI 14,130 km # https://webarchive.nationalarchives.gov.uk/20160108051201/http://www.ons.gov.uk/ons/guide-method/geography/beginner-s-guide/administrative/the-countries-of-the-uk/index.html

  df_long      <- melt(df, id=c("year"), variable.name = "land_use", value.name = "area")
  names(df_long) <- c("time", "u", "area" )
  notDuplictates <- (df_long$time - as.integer(df_long$time)) == 0
  df_long <- df_long[notDuplictates, ]
  df_long$country <- "N"
  df_A_NIr <- df_long
  return(list(df_A = df_A_NIr, df_A_orig = df_A_NIr))
}

## ---- combine_AgCensus

#' Function to wrangle AgCensus data 
#'  from text files to R objects
#'
#' @param l_df A list of data frames to combine
#' @return A BLAG object
#' @export
#' @examples
#' l_df = list(df_A_Eng, df_A_Sco, df_A_Wal, df_A_NIr)
#' x <- combine_AgCensus(l_df)
combine_AgCensus <- function(l_df = list(df_A_Eng, df_A_Sco, df_A_Wal, df_A_NIr)){
  df <- bind_rows(l_df)
  dt <- as.data.table(df)
  dt_A <- dt[, .(area = sum(area)), by = .(time, u)]
  dt_A$data_source <- "AgCensus"
  dt_A$area <- set_units(dt_A$area, km^2)
  dt_A$area <- set_units(dt_A$area, m^2)
  dt_A$area <- drop_units(dt_A$area)

  dt_D <- dt_A[, .(time, area = c(0, diff(area)), data_source), by = .(u)]
  dt_D <- dt_D[time >= 1900]
 
  return(list(dt_A = dt_A, dt_D = dt_D)) #, 
              # df_A_Eng = l_df$df_A_Eng, df_A_Sco = l_df$df_A_Sco, 
              # df_A_Wal = l_df$df_A_Wal, df_A_NIr = l_df$df_A_NIr))
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
  dt <- rbindlist(lapply(l_blags, '[[', "dt_B"), use.names=TRUE)
  dt$area <-  set_units(dt$area, m^2)
  dt$area <-  set_units(dt$area, km^2)
  dt$area <- drop_units(dt$area)
  dt_B <- dt
  
  # re-ordering by (u_to, u_from) (not u_from, u_to) allows 
  # the default vector to matrix conversion to work (where byrow = FALSE)
  dt_B <- arrange(dt_B, time, data_source, u_to, u_from)

  # G time series
  dt <- rbindlist(lapply(l_blags, '[[', "dt_G"), use.names=TRUE)
  dt$area <-  set_units(dt$area, m^2)
  dt$area <-  set_units(dt$area, km^2)
  dt$area <- drop_units(dt$area)
  dt_G <- dt
  
  # L time series
  dt <- rbindlist(lapply(l_blags, '[[', "dt_L"), use.names=TRUE)
  dt$area <-  set_units(dt$area, m^2)
  dt$area <-  set_units(dt$area, km^2)
  dt$area <- drop_units(dt$area)
  dt_L <- dt  
  
  # D time series
  dt <- rbindlist(lapply(l_blags, '[[', "dt_D"), use.names=TRUE)
  dt$area <-  set_units(dt$area, m^2)
  dt$area <-  set_units(dt$area, km^2)
  dt$area <- drop_units(dt$area)
  dt_D <- dt

  # A time series
  dt <- rbindlist(lapply(l_blags, '[[', "dt_A"), use.names=TRUE)
  dt$area <-  set_units(dt$area, m^2)
  dt$area <-  set_units(dt$area, km^2)
  dt$area <- drop_units(dt$area)
  dt_A <- dt
  
  return(list(dt_A = dt_A, dt_D = dt_D, dt_B = dt_B, dt_G = dt_G, dt_L = dt_L))
}


## ---- set_exclusions

#' Function to combine observations
#'  in BLAG objects to produce a single data table
#'
#' @param l_blags List of blag objects to combine
#' @return A blag object
#' @export
#' @examples
#' obs <- tar_read(c_obs)
#' unique(obs$dt_A$data_source[obs$dt_A$u == "other"])
#' obs <- set_exclusions(obs)
#' unique(obs$dt_A$data_source[obs$dt_A$u == "other"])
set_exclusions <- function(obs){

  # zeroes are actually missing values
  # should remove earlier
  obs$dt_A$area[obs$dt_A$area == 0]    <- NA
  obs$dt_D$area[obs$dt_D$area == 0]    <- NA
  obs$dt_D$area[is.nan(obs$dt_D$area)] <- NA
  obs$dt_G$area[obs$dt_G$area == 0]    <- NA
  obs$dt_G$area[is.nan(obs$dt_G$area)] <- NA
  obs$dt_L$area[obs$dt_L$area == 0]    <- NA
  obs$dt_L$area[is.nan(obs$dt_L$area)] <- NA
  obs$dt_B$area[obs$dt_B$area == 0] <- NA

  unique(obs$dt_D$data_source); unique(obs$dt_B$data_source)

  obs$dt_A$useData <- TRUE
  obs$dt_G$useData <- TRUE
  obs$dt_L$useData <- TRUE
  obs$dt_D$useData <- TRUE

  # exclude these data sources for woods
  data_sources_toExclude <- c("AgCensus", "IACS", "LCC", "CROME")
  obs$dt_A$useData[obs$dt_A$u == "woods" & obs$dt_A$data_source %in% data_sources_toExclude] <- FALSE
  obs$dt_G$useData[obs$dt_G$u == "woods" & obs$dt_G$data_source %in% data_sources_toExclude] <- FALSE
  obs$dt_L$useData[obs$dt_L$u == "woods" & obs$dt_L$data_source %in% data_sources_toExclude] <- FALSE
  obs$dt_D$useData[obs$dt_D$u == "woods" & obs$dt_D$data_source %in% data_sources_toExclude] <- FALSE

  # exclude these data sources for grass
  data_sources_toExclude <- c("IACS")
  obs$dt_A$useData[obs$dt_A$u == "grass" & obs$dt_A$data_source %in% data_sources_toExclude] <- FALSE

  # exclude these data sources for rough
  data_sources_toExclude <- c("IACS", "LCC")
  obs$dt_A$useData[obs$dt_A$u == "rough" & obs$dt_A$data_source %in% data_sources_toExclude] <- FALSE
    
  # exclude these data sources for urban
  obs$dt_A$useData[obs$dt_A$u == "urban" & obs$dt_A$data_source %in% data_sources_toExclude] <- FALSE
  obs$dt_G$useData[obs$dt_G$u == "urban" & obs$dt_G$data_source %in% data_sources_toExclude] <- FALSE
  obs$dt_L$useData[obs$dt_L$u == "urban" & obs$dt_L$data_source %in% data_sources_toExclude] <- FALSE
    
  # exclude these data sources for other
  obs$dt_A$useData[obs$dt_A$u == "other" & obs$dt_A$data_source %in% data_sources_toExclude] <- FALSE
  obs$dt_G$useData[obs$dt_G$u == "other" & obs$dt_G$data_source %in% data_sources_toExclude] <- FALSE
  obs$dt_L$useData[obs$dt_L$u == "other" & obs$dt_L$data_source %in% data_sources_toExclude] <- FALSE

  obs$dt_D <- subset(obs$dt_D, useData)
  obs$dt_A <- subset(obs$dt_A, useData)
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
    v_llik_B <- dnorm(dt_B$area, mean = dt_B$pred, sd = dt_B$sigma*abs(dt_B$area), log = T)
    v_llik_G <- dnorm(dt_G$area, mean = dt_G$pred, sd = dt_G$sigma*abs(dt_G$area), log = T)
    v_llik_L <- dnorm(dt_L$area, mean = dt_L$pred, sd = dt_L$sigma*abs(dt_L$area), log = T)
    v_llik_D <- dnorm(dt_D$area, mean = dt_D$pred, sd = dt_D$sigma*abs(dt_D$area), log = T)
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
#' x <- run_lcm_job(fname_job)
run_mcmc_job <- function(fname_job = "./slurm/run_mcmc_beta.job"){
  cmd <- paste0("sbatch ", fname_job)
  # submit the jobs to the SLURM queue
  err <- system(cmd)
  # and return the years and paths of the output files
  # these need to match slurm/process_LCM.R - no programmed check they are consistent
  v_times <- c(1990, 2015, 2017, 2018, 2019)
  return(list(
    v_times = v_times,
    v_fnames = paste0("./data/MCMC/")
  ))
}


## ---- get_rmse

#' Function to run MCMC processing job 
#'  for Beta matrix
#'
#' @param df A data frame containing the variables
#' @param v A character string for the name of the test variable
#' @param v_ref A character string for the name of the reference variable
#' @return Numeric The root-mean-square error
#' @export
#' @examples
#' rmse <- get_rmse(df = df, v = "IACS", v_ref = "Ref")
get_rmse <- function(df = df, v, v_ref = "Ref"){
  v     <- df[[v]]
  v_ref <- df[[v_ref]]
  resid <- v - v_ref
  rmse <- sqrt(mean(resid^2, na.rm = TRUE))
  # if no data, these will be NaN, which need to be NA
  rmse[is.nan(rmse)] <- NA
  return(rmse)
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
#' df <- get_uncert_scaling(obs, v_names_sources = c("AgCensus", "CS", "FC", "LCM", "CORINE", "LCC", "IACS", "CROME"))
get_uncert_scaling <- function(obs, v_names_sources = 
  c("AgCensus", "CS", "FC", "LCM", "CORINE", "LCC", "IACS", "CROME"),
  cv_AgCensus = 0.1){
  
  dt_D <- obs$dt_D
  df <- pivot_wider(dt_D, names_from = data_source, values_from = area)
  df <- subset(df, time >= 1990 & time <= 2020)
  df$Ref <- df$AgCensus
  df$Ref[df$u == "woods"] <- df$FC[df$u == "woods"]

  v_rmse <- sapply(v_names_sources, get_rmse, df = df, v_ref = "Ref")
  v_r2   <- sapply(v_names_sources, get_r2,   df = df, v_ref = "Ref")

  df <- data.frame(RMSE = v_rmse, r2 = v_r2, 
    # reduce RMSE proportional to r2, so that abs and prop measures contribute to sigma weighting
    sigma = v_rmse * abs(1 - v_r2))

  df <- df[order(df$sigma),]
  # AgCensus and FC form the reference, so are rows 1:2 when ordered
  # guess sigma for these as half the lowest value, which will be row 3 
  # very arbitrary assumption, to be improved upon
  df["AgCensus",]$sigma <- df[3,]$sigma * 0.5
  df["FC",]$sigma       <- df[3,]$sigma * 0.5

  df$sigma <- df$sigma / df["AgCensus",]$sigma * cv_AgCensus

  # add dummy values for false positive and neg rates
  df$Fp <- 0
  df$Fn <- 0
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
#' df <- add_uncert(obs, v_names_sources = c("AgCensus", "CS", "FC", "LCM", "CORINE", "LCC", "IACS", "CROME"))
add_uncert <- function(obs, df_uncert){
  
  df_uncert$data_source <- rownames(df_uncert)
  df_uncert <- df_uncert[, c("data_source", "sigma", "Fp", "Fn")]
  obs$dt_B <- merge(obs$dt_B, df_uncert, all.x = TRUE, by = "data_source") 
  obs$dt_G <- merge(obs$dt_G, df_uncert, all.x = TRUE, by = "data_source") 
  obs$dt_L <- merge(obs$dt_L, df_uncert, all.x = TRUE, by = "data_source") 
  obs$dt_A <- merge(obs$dt_A, df_uncert, all.x = TRUE, by = "data_source") 
  obs$dt_D <- merge(obs$dt_D, df_uncert, all.x = TRUE, by = "data_source") 
  
  return(obs)
}

