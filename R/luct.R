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

## ---- wrangle_CS

#' Function to wrangle Countryside Survey data 
#'  from text files to R objects
#'
#' @param fpath Filepath to text file
#' @return A BLAG object
#' @export
#' @examples
#' blag_CS <- wrangle_CS(fpath = "../data-raw/CS/UK_LUC_matrices_2018i.csv")
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
  return(list(dt_B_cs, dt_L_cs, dt_dA_cs, dt_G_cs))
}

# readCS()
readCS <- function(fpath = "../data-raw/CS/UK_LUC_matrices_2018i.csv"){
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

  p <- ggplot(data = subset(dt_B, year == 1991), aes(u_to, u_from)) 
  p <- p + geom_raster(aes(fill = area))
  p <- p + geom_text(aes(label = floor(area)), colour = "blue")
  p <- p + ylab("Land use in 1990")
  p <- p + xlab("Land use in 1991")
  p <- p + facet_wrap(~ year)
  p + scico::scale_fill_scico(palette = "lajolla")


  ## ----CSplotB, eval=TRUE, echo=FALSE, warning=FALSE, message=FALSE, fig.cap = "Time series of $\\mathbf{B}$ matrix, showing the areas changing land use over time. The layout of panels follows the matrix itself, so rows represent the starting land use, columns represent the end land use. We assumed that the rates of change were constant during the period between surveys."----
  p <- ggplot(dt_B, aes(year, area))
  p <- p + geom_line()
  p <- p + geom_point()
  p <- p + ylab(expression(paste(Area*", "*~km^2/y)))
  p <- p + facet_grid(u_from ~ u_to, scales = "fixed")
  p <- p + scale_x_continuous(limits = c(1970, 2020), breaks=c(1980, 2010))
  p


  ## ----CSplotG, eval=TRUE, echo=FALSE, warning=FALSE, message=FALSE, fig.cap = "Time series of implied area gains $\\mathbf{G}$ to each land use, from CS data. We assumed that the rates of change were constant during the period between surveys."----
  # calc gross gain from B; MARGIN = 3 because that is the time dimension
  dt_G <- as.data.table(t(apply(a_B[,,], MARGIN = 3, FUN = getAreaGrossGain_fromBeta, n_u = 6)))
  names(dt_G) <- names_u
  dt_G <- data.table(time = as.numeric(rownames(dt_G)), dt_G)
  dt_G <- melt(dt_G, id=c("time"), variable.name = "u", value.name = "area")
  dt_G$year <- v_times[dt_G$time]

  p <- ggplot(dt_G, aes(year, area))
  p <- p + geom_line()
  p <- p + geom_point()
  p <- p + ylab(expression(paste(Area*", "*~km^2/y)))
  p <- p + ggtitle("Gross Gains")
  p <- p + xlim(1970, 2020)
  p <- p + facet_wrap(~ u, scales = "fixed")
  p


  ## ----CSplotL, eval=TRUE, echo=FALSE, warning=FALSE, message=FALSE, fig.cap = "Time series of implied losses of area $\\mathbf{L}$ from each land use, from CS data. We assumed that the rates of change were constant during the period between surveys."----
  # calc gross loss from B; MARGIN = 3 because that is the time dimension
  dt_L <- as.data.table(t(apply(a_B[,,], MARGIN = 3, FUN = getAreaGrossLoss_fromBeta, n_u = 6)))
  names(dt_L) <- names_u
  dt_L <- data.table(time = as.numeric(rownames(dt_L)), dt_L)
  dt_L <- melt(dt_L, id=c("time"), variable.name = "u", value.name = "area")
  dt_L$year <- v_times[dt_L$time]

  p <- ggplot(dt_L, aes(year, area))
  p <- p + geom_line()
  p <- p + geom_point()
  p <- p + ggtitle("Gross Losses")
  p <- p + ylab(expression(paste(Area*", "*~km^2/y)))
  p <- p + xlim(1970, 2020)
  p <- p + facet_wrap(~ u, scales = "fixed")
  p


  ## ----CSplotdA, eval=TRUE, echo=FALSE, warning=FALSE, message=FALSE, fig.cap = "Time series of implied net change in area $\\Delta \\mathbf{A}$ of each land use, from CS data. We assumed that the rates of change were constant during the period between surveys."----
  # dt_dA_cs, ## this is dA not A**
  # calc net change from B; MARGIN = 3 because that is the time dimension
  dt_dA <- as.data.table(t(apply(a_B[,,], MARGIN = 3, FUN = getAreaNetChange_fromBeta, n_u = 6)))
  names(dt_dA) <- names_u
  dt_dA <- data.table(time = as.numeric(rownames(dt_dA)), dt_dA)
  dt_dA <- melt(dt_dA, id=c("time"), variable.name = "u", value.name = "area")
  dt_dA$year <- v_times[dt_dA$time]

  p <- ggplot(dt_dA, aes(year, area))
  p <- p + geom_line()
  p <- p + geom_point()
  p <- p + ylab(expression(paste(Area*", "*~km^2/y)))
  p <- p + ggtitle("Net Change")
  p <- p + xlim(1970, 2020)
  p <- p + facet_wrap(~ u, scales = "fixed")
  p


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
  return(list(dt_B_cs, dt_L_cs, dt_dA_cs, dt_G_cs))
}
