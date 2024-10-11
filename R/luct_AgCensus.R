v_region <- c("en", "sc", "wa", "ni", "uk")

  joinup <- function(df1, df2){
    # names of the columns to adjust
    colnames_to_adjust <- c("woods", "crops", "grass", "rough")

    # get the ending values of the first df
    v_area_end <- df1[length(df1$time), colnames_to_adjust]
    # get the starting values of the second df
    v_area_start <- df2[1, colnames_to_adjust]
    v_diff <- v_area_start - v_area_end

    # can't adjust NAs
    v_diff[is.na(v_diff)] <- 0

    # apply this difference to the first df
    df1$woods <- df1$woods + v_diff$woods
    df1$crops <- df1$crops + v_diff$crops
    df1$grass <- df1$grass + v_diff$grass
    df1$rough <- df1$rough + v_diff$rough

    return(df1)
}

#' @examples
#' agc_hist <- wrangle_AgCensus_historical(fname = "data-raw/AgCensus/all_LU_stats_FULLSERIES_km2.csv")
#' str(agc_hist)
wrangle_AgCensus_historical <- function(
  fname = "data-raw/AgCensus/all_LU_stats_FULLSERIES_km2.csv"){

  df <- read.csv(fname)

  # for consistency, just use first two characters, incl capitalisation
  df$region <- tolower(substr(df$country, 1, 2))
  df$country <- NULL

  # sum to UK in a separate data frame
  # remove "region" before summing all - otherwise causes an error
  df_uk <- subset(df, select = -region) # WILL
  df_uk <- df_uk %>% group_by(year) %>%
      summarise_all(list(~ sum(., na.rm=TRUE)))
  df_uk$region <- "uk"

  # add this back in to the original
  df <- rbind(df_uk, df)
  df$region <- factor(df$region, levels = v_region)

  # and declare units as km2
  df$crops <- set_units(df$crops, km^2)
  df$grass <- set_units(df$grass, km^2)
  df$rough <- set_units(df$rough, km^2)
  df$woods <- set_units(df$woods, km^2)

  # make names and their order consistent
  names(df)[names(df) == "year"] <- "time"
  names(df)[names(df) == "land_use"] <- "u"
  df <- df[, c("region", "time", "woods", "crops", "grass", "rough" )]

  # and make a long-format version
  df_wide <- df
  df_long      <- melt(df_wide, id=c("time", "region"), variable.name = "u", value.name = "area")

  return(list(df_long = df_long, df_wide = df_wide))
}


## ---- wrangle_AgCensus_Eng_1

#' Function to wrangle AgCensus data
#'  from text files to R objects
#'
#' @param v_fpath Vector of filepaths to data files
#' @return A BLAG object
#' @export
#' @examples
#' fpath = "./data-raw/AgCensus/England/AgCensus_England_ha_1900-2010.csv"
#' df1 <- wrangle_AgCensus_Eng_1(fpath1)
wrangle_AgCensus_Eng_1 <- function(fpath =
  "./data-raw/AgCensus/England/AgCensus_England_ha_1900-2010.csv"){

#<!--- { read_data_England1 -->
  df <- read.csv(fpath, na.strings = c("n/c", "-"))
  names(df)[names(df) == "year"] <- "time"

  df$woods <-NA
  df$crops <- # crops are in columns 3-28
    rowSums(df[, 3:28], na.rm = TRUE)
  df$grass <-
    df$Permanent.grass..excl.rough.grazing.
  df$rough <-
    df$Total.rough.grazing
  df <- df[, c("time", "woods", "crops", "grass", "rough" )]

  # interpolate between decadal points
  df_allyears <- data.frame(time = 1900:1990)
  df <- merge(df_allyears, df, all.x = TRUE)
  df[,-1] <- na.approx(df[,-1])

  # set units to hectares (displayed in Excel as kha, but stored as ha)
  df[,-1] <- df[,-1] %>% map2_dfc("ha", ~set_units(.x, .y, mode = "standard"))
  # and convert units to km2
  df[,-1] <- df[,-1] %>% map2_dfc("km^2", ~set_units(.x, .y, mode = "standard"))

  df_hist <- df # make a copy of original data
  return(df)
}


## ---- wrangle_AgCensus_Eng_2

#' Function to wrangle AgCensus data
#'  from text files to R objects
#'
#' @param v_fpath Vector of filepaths to data files
#' @return A BLAG object
#' @export
#' @examples
#' fpath = "./data-raw/AgCensus/England/AgCensus_England_ha_1983-2019.csv"
#' df2 <- wrangle_AgCensus_Eng_2(fpath)
wrangle_AgCensus_Eng_2 <- function(fpath =
  "./data-raw/AgCensus/England/AgCensus_England_ha_1983-2019.csv"){

#<!--- { read_data_England2 -->
  df <- read.csv(fpath, na.strings = c("n/c"))
  names(df)[names(df) == "year"] <- "time"

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
  df <- df[, c("time", "woods", "crops", "grass", "rough" )]

  # code the first time as a known change point (better than in the data file), but still ugly
  df$time[1] <- df$time[1] + 0.5

  # set units to hectares (displayed in Excel as kha, but stored as ha)
  df[,-1] <- df[,-1] %>% map2_dfc("ha", ~set_units(.x, .y, mode = "standard"))
  # and convert units to km2
  df[,-1] <- df[,-1] %>% map2_dfc("km^2", ~set_units(.x, .y, mode = "standard"))

  # remove overlap between hist and contemporary data, which starts in 1983
  # remove overlap with contemporary data
  # df_hist <- subset(df_hist, time <= 1983)
  # df <- rbind(df_hist, df)

  return(df)
}


## ---- wrangle_AgCensus_Sco_1

#' Function to wrangle AgCensus data
#'  from text files to R objects
#'
#' @param v_fpath Vector of filepaths to data files
#' @return A BLAG object
#' @export
#' @examples
#' fpath = "./data-raw/AgCensus/Scotland/AgCensus_Scotland_ha_1883-2014.csv"
#' df1 <- wrangle_AgCensus_Sco_1(fpath)
wrangle_AgCensus_Sco_1 <- function(fpath =
  "./data-raw/AgCensus/Scotland/AgCensus_Scotland_ha_1883-2014.csv"){

  # read historical data
  df <- read.csv(fpath, na.strings = c("n/a", ":", "b"))
  names(df)[names(df) == "year"] <- "time"

  # set units to hectares (displayed in Excel as kha, but stored as ha)
  df[,-1] <- df[,-1] %>% map2_dfc("ha",   ~set_units(.x, .y, mode = "standard"))
  # and convert units to km2
  df[,-1] <- df[,-1] %>% map2_dfc("km^2", ~set_units(.x, .y, mode = "standard"))

  return(df)
}


## ---- wrangle_AgCensus_Sco_2

#' Function to wrangle AgCensus data
#'  from text files to R objects
#'
#' @param v_fpath Vector of filepaths to data files
#' @return A BLAG object
#' @export
#' @examples
#' fpath = "./data-raw/AgCensus/Scotland/AgCensus_Scotland_ha_2009-2020.csv"
#' df2 <- wrangle_AgCensus_Sco_2(fpath)
wrangle_AgCensus_Sco_2 <- function(fpath =
    "./data-raw/AgCensus/Scotland/AgCensus_Scotland_ha_2009-2019.csv"){

  # read current data
  df <- read.csv(fpath, na.strings = c("n/a", ":", "b"))
  names(df)[names(df) == "year"] <- "time"

  df$woods <-
    df$Other.land..including.woodland..3.
  df$crops <-
    df$Total.crops..fallow..and.set.aside.
  df$grass <-
    df$Total.grass
  df$rough <-
    df$Rough.grazing
  df <- df[, c("time", "woods", "crops", "grass", "rough" )]

  # set units to hectares (displayed in Excel as kha, but stored as ha)
  df[,-1] <- df[,-1] %>% map2_dfc("ha",   ~set_units(.x, .y, mode = "standard"))
  # and convert units to km2
  df[,-1] <- df[,-1] %>% map2_dfc("km^2", ~set_units(.x, .y, mode = "standard"))

  return(df)
}


## ---- wrangle_AgCensus_Wal_1

#' Function to wrangle AgCensus data
#'  from text files to R objects
#'
#' @param v_fpath Vector of filepaths to data files
#' @return A BLAG object
#' @export
#' @examples
#'  fpath1 = "./data-raw/AgCensus/Wales/AgCensus_Wales_ha_1867-2012.csv"
#'  fpath2 = "./data-raw/AgCensus/Wales/AgCensus_Wales_ha_1998-2019.csv"
#'  x <- wrangle_AgCensus_Wal_1(fpath)
wrangle_AgCensus_Wal_1 <- function(fpath =
  "./data-raw/AgCensus/Wales/AgCensus_Wales_ha_1867-2012.csv"){

  # read historical data
  df <- read.csv(fpath, na.strings = c("n/a", ":", "b"))
  names(df)[names(df) == "year"] <- "time"

  df$woods <- NA
  df$crops <-
    df$Crops
  df$grass <-
    df$Permanent.pasture +
    df$New.grass
  df$rough <-
    df$Rough.grazing
  df <- df[, c("time", "woods", "crops", "grass", "rough" )]

  # set units to hectares (displayed in Excel as kha, but stored as ha)
  df[,-1] <- df[,-1] %>% map2_dfc("ha",   ~set_units(.x, .y, mode = "standard"))
  # and convert units to km2
  df[,-1] <- df[,-1] %>% map2_dfc("km^2", ~set_units(.x, .y, mode = "standard"))

  return(df)
}


## ---- wrangle_AgCensus_Wal_2

#' Function to wrangle AgCensus data
#'  from text files to R objects
#'
#' @param v_fpath Vector of filepaths to data files
#' @return A BLAG object
#' @export
#' @examples
#'  fpath = "./data-raw/AgCensus/Wales/AgCensus_Wales_ha_1998-2019.csv"
#'  x <- wrangle_AgCensus_Wal_2(fpath)
wrangle_AgCensus_Wal_2 <- function(fpath =
    "./data-raw/AgCensus/Wales/AgCensus_Wales_ha_1998-2019.csv"){

  # read current data
  df <- read.csv(fpath, na.strings = c("n/a"))
  names(df)[names(df) == "year"] <- "time"

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
  df <- df[, c("time", "woods", "crops", "grass", "rough" )]

  # set units to hectares (displayed in Excel as kha, but stored as ha)
  df[,-1] <- df[,-1] %>% map2_dfc("ha",   ~set_units(.x, .y, mode = "standard"))
  # and convert units to km2
  df[,-1] <- df[,-1] %>% map2_dfc("km^2", ~set_units(.x, .y, mode = "standard"))

  return(df)
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
  fpath = "./data-raw/AgCensus/NIreland/AgCensus_NIreland_ha_1981-2019.csv",
  agc_hist){

  # read NIr data
  df <- read.csv(fpath, na.strings = c("n/a", ":", "b"))
  names(df)[names(df) == "year"] <- "time"

  df$woods <- df$Woods.and.plantations
  df$crops <-
    df$TOTAL.CROPS +
    df$Set.aside
  df$grass <- df$TOTAL.GRASS
  df$rough <- df$Hill.or.Rough.Land
  df <- df[, c("time", "woods", "crops", "grass", "rough" )]

  # set units to hectares
  df[,-1] <- df[,-1] %>% map2_dfc("ha",   ~set_units(.x, .y, mode = "standard"))
  # and convert units to km2
  df[,-1] <- df[,-1] %>% map2_dfc("km^2", ~set_units(.x, .y, mode = "standard"))
  return(df)
}

## ---- join_AgCensus

#' Function to concatenate two or three AgCensus data
#'  frames together, removing overlap and
#'  correcting for the difference at the join
#'
#' @param df0 Data frame for first period
#' @param df1 Data frame for first period
#' @param df2 Data frame for first period
#' @return A concatenated data frame for first period
#' @export
#' @examples
#' df <- join_AgCensus(agc_hist$df_wide, df1, df2)
join_AgCensus <- function(df0, df1, df2){

  df0 <- df0[, c("time", "woods", "crops", "grass", "rough" )]
  df1 <- df1[, c("time", "woods", "crops", "grass", "rough" )]
  df2 <- df2[, c("time", "woods", "crops", "grass", "rough" )]
  # set the order
  df0 <- dplyr::arrange(df0, time)
  df1 <- dplyr::arrange(df1, time)
  df2 <- dplyr::arrange(df2, time)

  summary(df0$time)
  summary(df1$time)
  summary(df2$time)

  # remove overlap, keeping the more recent data
  start_df2 <- df2$time[1]
  df0 <- subset(df0, time <= start_df2)
  df1 <- subset(df1, time <= start_df2)

  start_df1 <- df1$time[1]
  df0 <- subset(df0, time <= start_df1)

  df1 <- joinup(df1, df2)
  df0 <- joinup(df0, df1)

  df <- rbind(df0, df1, df2)

  return(df)
}

#' function not used
# dfs <- remove_steps(df,
  # knownChangePoint = NULL,
  # unknownChangePoint = c(1900))
remove_steps <- function(df,
  knownChangePoint = c(1983, 2009),
  ## 1977 and 2006 appear to be step changes but there is no record of a change in method
  unknownChangePoint = c(1900) # NULL or add apparent step changes here
  ){

  # need to repeat for all variables
  df <- drop_units(df) # gam doesn't work with units
  for (j in 2:(length(df))){   # for each variable except year
  #j = 2
    # go through all the columns except first (= time)
    # first, adjust for step changes where values are given under old and new basis
    for (i in length(knownChangePoint):1){
    #i = 1
      df_pre  <- subset(df, time <= knownChangePoint[i])
      df_post <- subset(df, time > knownChangePoint[i])
      obs_pre  <- df[,j][df$time == knownChangePoint[i]]
      obs_post <- df[,j][df$time == (knownChangePoint[i]+0.5)]
      jump <- obs_post - obs_pre
      df_pre[,j] <- df_pre[,j] + jump
      df <- rbind(df_pre, df_post)
    } # for each knownChangePoint

    # gam version for predicting area in changepoint year, based on a gam fitted to the previous data
    # use this where there are step changes, but values not given on old and new basis
    if (length(unknownChangePoint) > 0){
      for (i in length(unknownChangePoint):1){
      #i = 2
        df_pre  <- subset(df, time < unknownChangePoint[i])
        df_post <- subset(df, time >= unknownChangePoint[i])
        if ( sum(!is.na(df_pre[,j])) > 10 ){ # only do if >10 values, excluding NAs
          #k <- max(length(df_pre[,1])-1, 4)
          k <- min(length(df_pre[,1])-1, 10)
          mod_pre  <- gam(df_pre[,j] ~ s(time, bs = "cr", k = k), data = df_pre)
          #k <- max(length(df_post[,1])-1, 4)
          #k <- min(length(df_post[,1])-1, 10)
          #mod_post <- gam(grassArea ~ s(time, bs = "cr", k = k), data = df_post)
          pred_cpt_pre  <- predict.gam(mod_pre, data.frame(time=unknownChangePoint[i]))
          #pred_cpt_post <- predict.gam(mod_post, data.frame(time=unknownChangePoint[i]))
          #pred_pre  <- predict.gam(mod_pre, df)
          obs_pre  <- df[,j][df$time == unknownChangePoint[i]]
          jump <- obs_pre - pred_cpt_pre
          df_pre[,j] <- df_pre[,j] + jump
          df <- rbind(df_pre, df_post)
        } # >10 values
      }   # i change points
    }     # if there are any change points
  }       # for each variable except time

  # remove duplicated rows with knownChangePoints, where one has time = XXXX.5
  #df <- subset(df, time %% floor(time) == 0) # so modulo division is not zero
  return(df)
}


process_AgCensus_Eng <- function(
  fpath = "./data-raw/AgCensus/England/AgCensus_England_ha_1983-2019.csv",
  agc_hist
  ) {

  df0 <- subset(agc_hist$df_wide, region == "en")
  # df1 <- wrangle_AgCensus_Eng_1(v_fpath[1])
  df2 <- wrangle_AgCensus_Eng_2(fpath)

  # fill in the missing woods
  # get the ending values of the first df
  v_area_end <- df0$woods[length(df0$time)]
  # get the starting values of the second df
  v_area_start <- df2$woods[1]
  v_diff <- v_area_start - v_area_end
  # apply this difference to the second df
  df2$woods <- df2$woods - v_diff

  df <- join_AgCensus(df0, df2, df2)


  # estimate the nadir of woods area
  df$woods[df$time == 1940] <- 6000
  # and fill in the gaps in woods
  df$woods <- na.approx(df$woods)
  df$rough <- na.approx(df$rough)

  # make a long format version
  df_long      <- melt(df, id=c("time"), variable.name = "u", value.name = "area")
  notDuplictates <- (df_long$time - as.integer(df_long$time)) == 0
  df_long <- df_long[notDuplictates, ]
  df_long$region <- "en"
  df$region      <- "en"

  return(list(df_long = df_long, df_wide = df))
}

process_AgCensus_Sco <- function(v_fpath = c(
  "./data-raw/AgCensus/Scotland/AgCensus_Scotland_ha_1883-2014.csv",
  "./data-raw/AgCensus/Scotland/AgCensus_Scotland_ha_2009-2019.csv"),
  agc_hist
  ){

  df0 <- subset(agc_hist$df_wide, region == "sc")
  df0 <- df0[, c("time", "woods", "crops", "grass", "rough" )]
  df1 <- wrangle_AgCensus_Sco_1(v_fpath[1])
  df2 <- wrangle_AgCensus_Sco_2(v_fpath[2])

  df <- join_AgCensus(df0, df1, df2)

  # fill in the missing rough values around 1960
  df$rough[df$time == 1959] <- NA
  df$rough[df$time == 1960] <- NA

  # estimate the nadir of woods area
  df$woods[df$time == 1940] <- 1000
  # and fill in the gaps in woods
  df$woods <- na.approx(df$woods)
  df$rough <- na.approx(df$rough)

  # make a long format version
  df_long      <- melt(df, id=c("time"), variable.name = "u", value.name = "area")
  notDuplictates <- (df_long$time - as.integer(df_long$time)) == 0
  df_long <- df_long[notDuplictates, ]
  df_long$region <- "sc"
  df$region      <- "sc"

  return(list(df_long = df_long, df_wide = df))
}

process_AgCensus_Wal <- function(v_fpath = c(
 "./data-raw/AgCensus/Wales/AgCensus_Wales_ha_1867-2012.csv",
 "./data-raw/AgCensus/Wales/AgCensus_Wales_ha_1998-2019.csv"),
  agc_hist
  ){

  df0 <- subset(agc_hist$df_wide, region == "wa")
  df0 <- df0[, c("time", "woods", "crops", "grass", "rough" )]
  df1 <- wrangle_AgCensus_Wal_1(v_fpath[1])
  df2 <- wrangle_AgCensus_Wal_2(v_fpath[2])

  df <- join_AgCensus(df0, df1, df2)

  # estimate the nadir of woods area
  df$woods[df$time == 1940] <- 250
  # and fill in the gaps in woods
  df$woods <- na.approx(df$woods)
  df$rough <- na.approx(df$rough)

  # make a long format version
  df_long      <- melt(df, id=c("time"), variable.name = "u", value.name = "area")
  notDuplictates <- (df_long$time - as.integer(df_long$time)) == 0
  df_long <- df_long[notDuplictates, ]
  df_long$region <- "wa"
  df$region      <- "wa"

  return(list(df_long = df_long, df_wide = df))
}

process_AgCensus_NIr <- function(fpath =
 "./data-raw/AgCensus/Wales/AgCensus_Wales_ha_1867-2012.csv",
  agc_hist
  ){

  df0 <- subset(agc_hist$df_wide, region == "ni")
  df0 <- df0[, c("time", "woods", "crops", "grass", "rough" )]
  df1 <- wrangle_AgCensus_NIr(fpath)

  df <- join_AgCensus(df0, df1, df1)

  # interpolate between the gap from 1969 to 1983
  df_allyears <- data.frame(time = 1750:2020)
  df <- merge(df_allyears, df, all.x = TRUE)
  df[,-1] <- na.approx(df[,-1])

  # make a long format version
  df_long      <- melt(df, id=c("time"), variable.name = "u", value.name = "area")
  notDuplictates <- (df_long$time - as.integer(df_long$time)) == 0
  df_long <- df_long[notDuplictates, ]
  df_long$region <- "ni"
  df$region      <- "ni"

  return(list(df_long = df_long, df_wide = df))
}

#' l_df = list(df_A_Eng, df_A_Sco, df_A_Wal, df_A_NIr)
#' x <- combine_AgCensus(l_df, agc_hist$df_long)
combine_AgCensus <- function(
  l_df = list(df_A_Eng, df_A_Sco, df_A_Wal, df_A_NIr),
  v_region = "uk"){

  #df <- bind_rows(l_df)
  df <- plyr::rbind.fill(l_df)
  dt <- as.data.table(df)
  dt_A <- dt[, .(area = mean(area, na.rm = TRUE)), by = .(region, time, u)]

  with(dt_A, table(region, time, u))
  dt_A_uk <- dt_A[region != "uk", .(area =  sum(area, na.rm = TRUE)), by = .(time, u)]
  dt_A_uk$region <- "uk"
  dt_A <- bind_rows(dt_A, dt_A_uk)
  dt_A$data_source <- "AgCensus"
  # we need to have all in m2 for consistency
  dt_A$area <- set_units(dt_A$area, km^2)
  dt_A$area <- set_units(dt_A$area, m^2)
  dt_A$area <- drop_units(dt_A$area)

  dt_D <- dt_A[, .(time, area = c(0, diff(area)), data_source), by = .(region, u)]
  dt_D <- dt_D[time >= 1900]

  # subset to return only the selected countries
  dt_A <- dt_A[region %in% v_region]
  dt_D <- dt_D[region %in% v_region]

  return(list(dt_A = dt_A, dt_D = dt_D)) #,
              # df_A_Eng = l_df$df_A_Eng, df_A_Sco = l_df$df_A_Sco,
              # df_A_Wal = l_df$df_A_Wal, df_A_NIr = l_df$df_A_NIr))
}


wrangle_MODIS_urban <- function(fpath =
  "./data-raw/MODIS/FAOStat_Artificial_land_surface_MODIS.csv"){

  # get the urban split by DA region in 2018
  fname <- here("data-raw/mask", "dt_land_100m.qs")
  dt_mask <- qread(fname)
  # str(dt_mask)
  # table(dt_mask$country)

  fname <- here("data/LCM/Level1", "r_U_lcm_100m_2019.tif")
  r <- raster(fname)
  # plot(r)
  dt <- as.data.table(as.data.frame(r))
  dt <- data.table(dt_mask, u = dt$r_U_lcm_100m_2019)
  dt_urban <- dt[u == 5, .N, by = country]
  dt_urban <- dt_urban[!is.na(country)]
  dt_urban[, frac_urban := N / sum(N)]
  dt_urban[, region := v_region[country]]

  # read in the MODIS data
  dt <- fread(fpath)
  dt_A <- with(dt, data.table(region = "uk", time = Year, area = Value, u = "urban", data_source = "MODIS"))
  dt_A$area <- set_units(dt_A$area*1000, ha)
  dt_A$area <- set_units(dt_A$area, m^2)
  dt_A$area <- drop_units(dt_A$area)

  dt_A_uk <- dt_A
  dt_A_en <- dt_A_uk
  dt_A_sc <- dt_A_uk
  dt_A_wa <- dt_A_uk
  dt_A_ni <- dt_A_uk

  dt_A_en$region <- "en"
  dt_A_sc$region <- "sc"
  dt_A_wa$region <- "wa"
  dt_A_ni$region <- "ni"

  dt_A_en[, area := area * dt_urban[region == "en", frac_urban]]
  dt_A_sc[, area := area * dt_urban[region == "sc", frac_urban]]
  dt_A_wa[, area := area * dt_urban[region == "wa", frac_urban]]
  dt_A_ni[, area := area * dt_urban[region == "ni", frac_urban]]

  dt_A <- rbindlist(list(
    dt_A_uk,
    dt_A_en,
    dt_A_sc,
    dt_A_wa,
    dt_A_ni
  ))

  # calculate net change  by difference
  dt_D <- dt_A[, .(time, area = c(NA, diff(area)), data_source), by = .(region, u)]
  # infer that change is all gain, no losses
  dt_G <- dt_D
  dt_L <- dt_D
  dt_L$area <- 0
  dim(dt_D)
  return(list(dt_G = dt_G, dt_L = dt_L, dt_A = dt_A, dt_D = dt_D))
}

wrangle_popcen <- function(fname = fs::path_rel(here("data-raw/pre1950/POPCEN",
      "dt_A_DevA_km2.qs"))) {

  dt <- qread(fname)

  setnames(dt, "country", "region")
  setnames(dt, "Year", "time")
  dt <- dt[, .(time, region, urban)]
  setnames(dt, "urban", "area")
  dt[, u := "urban"]
  dt[, data_source := "popcen"]
  dt <- dt[, .(region, u, time, area, data_source)]

  dt[region == "E", region := "en"]
  dt[region == "S", region := "sc"]
  dt[region == "W", region := "wa"]
  dt[region == "NI", region := "ni"]

  # set units
  units(dt$area)
  dt$area <- drop_units(dt$area)
  dt$area <- set_units(dt$area, km^2)
  dt$area <- set_units(dt$area, m^2)
  dt$area <- drop_units(dt$area)

  # fit a smooth curve to join popcen and lcm data without losses
  dt[, region := factor(region)]
  levels(dt$region)
  m <- gam(area ~ region + s(time, by = region, k = 4), data = dt)
  dt[, area_pred := predict(m)]
  # p <- ggplot(dt_G, aes(time, area, colour = region))
  # p <- p + geom_point()
  # p <- p + geom_line(aes(y= area_pred))
  # p <- p + facet_wrap(~ region, scale = "free_y")
  # p
  dt[, area := area_pred]
  dt[, area_pred := NULL]

  dt_A <- dt
# calculate net change  by difference
  dt_D <- dt_A[, .(time, area = c(NA, diff(area)), data_source), by = .(region, u)]
  # infer that change is all gain, no losses
  dt_G <- copy(dt_D)
  dt_L <- copy(dt_D)
  dt_L[, area := 0]
  dt_D[area < 0, area := 0]
  dt_G[area < 0, area := 0]
  return(list(dt_G = dt_G, dt_L = dt_L, dt_A = dt_A, dt_D = dt_D))
}

# tar_read(c_blag_MODIS, store = "_targets_wrangler")
