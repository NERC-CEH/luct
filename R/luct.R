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
#' dA_u <- getAreaNetChange_fromBeta(v_B, n_u = 6)
getAreaNetChange_fromBeta <- function(v_B, n_u = sqrt(length(v_B))){
  m_B <- matrix(v_B, n_u, n_u)
  diag(m_B) <- 0
  dA_u <- colSums(m_B) - rowSums(m_B)
  return(dA_u)  
}

## ---- wrangle_AgCensus_Eng

#' Function to wrangle AgCensus data 
#'  from text files to R objects
#'
#' @param fpath Filepath to text file
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
#' @param fpath Filepath to text file
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
#' @param fpath Filepath to text file
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
#' @param fpath Filepath to text file
#' @return A BLAG object
#' @export
#' @examples
#' l_df = list(df_A_Eng, df_A_Sco, df_A_Wal, df_A_NIr)
#' x <- combine_AgCensus(l_df)
combine_AgCensus <- function(l_df = list(df_A_Eng, df_A_Sco, df_A_Wal, df_A_NIr)){
  df <- bind_rows(l_df)
  dt <- as.data.table(df)
  dt_A_UK <- dt[, .(area = sum(area)), by = .(time, u)]
  dt_A_UK$data_source <- "AgCensus"

  dt_dA_UK <- dt_A_UK[, .(time, area = c(0, diff(area)), data_source), by = .(u)]
  dt_dA_UK <- dt_dA_UK[time >= 1900]
  dt_dA_UK$area <- set_units(dt_dA_UK$area, km^2)

  return(list(dt_A_UK = dt_A_UK, dt_dA_UK = dt_dA_UK)) #, 
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

  ## ----CSplotdA, eval=TRUE, echo=FALSE, warning=FALSE, message=FALSE, fig.cap = "Time series of implied net change in area $\\Delta \\mathbf{A}$ of each land use, from CS data. We assumed that the rates of change were constant during the period between surveys."----
  # dt_dA_cs, ## this is dA not A**
  # calc net change from B; MARGIN = 3 because that is the time dimension
  dt_dA <- as.data.table(t(apply(a_B[,,], MARGIN = 3, FUN = getAreaNetChange_fromBeta, n_u = 6)))
  names(dt_dA) <- names_u
  dt_dA <- data.table(time = as.numeric(rownames(dt_dA)), dt_dA)
  dt_dA <- melt(dt_dA, id=c("time"), variable.name = "u", value.name = "area")
  dt_dA$year <- v_times[dt_dA$time]

  ## ----saving, eval=TRUE, echo=FALSE, warning=FALSE, message=FALSE----------------------
  dt_dA$data_source <- "CS"
  dt_B$data_source <- "CS"
  dt_G$data_source <- "CS"
  dt_L$data_source <- "CS"

  dt_dA$time <- dt_dA$year; dt_dA$year <- NULL 
  dt_B$time <- dt_B$year; dt_B$year <- NULL 
  dt_G$time <- dt_G$year; dt_G$year <- NULL 
  dt_L$time <- dt_L$year; dt_L$year <- NULL

  dt_dA$area <- set_units(dt_dA$area, km^2)
  dt_B$area <- set_units(dt_B$area, km^2)
  dt_G$area <- set_units(dt_G$area, km^2)
  dt_L$area <- set_units(dt_L$area, km^2)

  dt_dA$area <- set_units(dt_dA$area, m^2)
  dt_B$area <- set_units(dt_B$area, m^2)
  dt_G$area <- set_units(dt_G$area, m^2)
  dt_L$area <- set_units(dt_L$area, m^2) 

  dt_dA$area <- drop_units(dt_dA$area)
  dt_B$area <- drop_units(dt_B$area)
  dt_G$area <- drop_units(dt_G$area)
  dt_L$area <- drop_units(dt_L$area) 

  dt_dA_cs <- dt_dA
  dt_B_cs <- dt_B
  dt_G_cs <- dt_G
  dt_L_cs <- dt_L 
  return(list(dt_B = dt_B, dt_L = dt_L, dt_dA = dt_dA, dt_G = dt_G))
}
