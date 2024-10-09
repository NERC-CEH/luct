library(targets)
library(tarchetypes) # for tar_knitr_deps_expr()
library(here) # construct file paths relative to project root
library(fs) # file system operations
library(future)
future::plan("multicore")

# All the functions to calculate the targets are defined here.
source(here::here("R", "luc_track.R"))
source(here::here("R", "luct.R"))
source(here::here("R", "luct_AgCensus.R"))

# Set target-specific options such as packages.
# It is recommended to load here all the packages required by functions managed only by {targets}.
options(tidyverse.quiet = TRUE)
options(bitmapType='cairo')
tar_option_set(
  packages = c(
    "dplyr", "purrr", "units", "data.table", "ggplot2", "stringr", "qs",
    "zoo", "mgcv", "reshape2", "readxl", "tidyr", "sp", "sf", "raster", #"rgeos", "ggforce", "plyr",
    "rgdal", "grid", "spCEH", "scico", "stars", "BayesianTools", "scales"),
  format = "qs",
  cue = tar_cue(
      mode = "thorough", #c("thorough", "always", "never"),
      command   = TRUE,
      depend    = TRUE,
      format    = TRUE,
      iteration = TRUE,
      file      = TRUE)
)

v_times <- 1900:2020
res <- 10000

c_post_B_en <- tar_read(c_post_B_en, store = "_targets_plotB")
c_post_B_sc <- tar_read(c_post_B_sc, store = "_targets_plotB")
c_post_B_wa <- tar_read(c_post_B_wa, store = "_targets_plotB")
c_post_B_ni <- tar_read(c_post_B_ni, store = "_targets_plotB")

c_fname_dt_Umap_en <- tar_read(c_fname_dt_Umap_en, store = "_targets_postU")
c_fname_dt_Umap_sc <- tar_read(c_fname_dt_Umap_sc, store = "_targets_postU")
c_fname_dt_Umap_wa <- tar_read(c_fname_dt_Umap_wa, store = "_targets_postU")
c_fname_dt_Umap_ni <- tar_read(c_fname_dt_Umap_ni, store = "_targets_postU")

# list of target objects
list(

  # Get the en B matrix for mineral soil only
  tar_target(
    c_a_B_en_min,
    get_B_for_soil_type(
      dt = qread(c_fname_dt_Umap_en),
      res = res,
      my_region = "en",
      my_soil_type = 0,
      v_times = v_times ),
     cue = tar_cue(mode = "thorough")
  ),

  # Get the en B matrix for organic soil only
  tar_target(
    c_a_B_en_org,
    get_B_for_soil_type(
      dt = qread(c_fname_dt_Umap_en),
      res = res,
      my_region = "en",
      my_soil_type = 1,
      v_times = v_times),
     cue = tar_cue(mode = "thorough")
  ),

  # Get the sc B matrix for mineral soil only
  tar_target(
    c_a_B_sc_min,
    get_B_for_soil_type(
      dt = qread(c_fname_dt_Umap_sc),
      res = res,
      my_region = "sc",
      my_soil_type = 0,
      v_times = v_times),
     cue = tar_cue(mode = "thorough")
  ),

  # Get the sc B matrix for organic soil only
  tar_target(
    c_a_B_sc_org,
    get_B_for_soil_type(
      dt = qread(c_fname_dt_Umap_sc),
      res = res,
      my_region = "sc",
      my_soil_type = 1,
      v_times = v_times),
     cue = tar_cue(mode = "thorough")
  ),

  # Get the wa B matrix for mineral soil only
  tar_target(
    c_a_B_wa_min,
    get_B_for_soil_type(
      dt = qread(c_fname_dt_Umap_wa),
      res = res,
      my_region = "wa",
      my_soil_type = 0,
      v_times = v_times),
     cue = tar_cue(mode = "thorough")
  ),

  # Get the wa B matrix for organic soil only
  tar_target(
    c_a_B_wa_org,
    get_B_for_soil_type(
      dt = qread(c_fname_dt_Umap_wa),
      res = res,
      my_region = "wa",
      my_soil_type = 1,
      v_times = v_times),
     cue = tar_cue(mode = "thorough")
  ),

  # Get the ni B matrix for mineral soil only
  tar_target(
    c_a_B_ni_min,
    get_B_for_soil_type(
      dt = qread(c_fname_dt_Umap_ni),
      res = res,
      my_region = "ni",
      my_soil_type = 0,
      v_times = v_times),
     cue = tar_cue(mode = "thorough")
  ),

  # Get the ni B matrix for organic soil only
  tar_target(
    c_a_B_ni_org,
    get_B_for_soil_type(
      dt = qread(c_fname_dt_Umap_ni),
      res = res,
      my_region = "ni",
      my_soil_type = 1,
      v_times = v_times),
     cue = tar_cue(mode = "thorough")
  ),

  # Get the en uncertainties in B from df to array of matrices format
  tar_target(
    c_l_a_B_uncert_en,
    get_B_uncertainty_as_array(
      df = c_post_B_en$df_B,
      my_region = "en",
      v_times = unique(c_post_B_en$df_D$time)),
    cue = tar_cue(mode = "thorough")
  ),

  # Get the sc uncertainties in B from df to array of matrices format
  tar_target(
    c_l_a_B_uncert_sc,
    get_B_uncertainty_as_array(
      df = c_post_B_sc$df_B,
      my_region = "sc",
      v_times = unique(c_post_B_sc$df_D$time)),
    cue = tar_cue(mode = "thorough")
  ),

  # Get the wa uncertainties in B from df to array of matrices format
  tar_target(
    c_l_a_B_uncert_wa,
    get_B_uncertainty_as_array(
      df = c_post_B_wa$df_B,
      my_region = "wa",
      v_times = unique(c_post_B_wa$df_D$time)),
    cue = tar_cue(mode = "thorough")
  ),

  # Get the ni uncertainties in B from df to array of matrices format
  tar_target(
    c_l_a_B_uncert_ni,
    get_B_uncertainty_as_array(
      df = c_post_B_ni$df_B,
      my_region = "ni",
      v_times = unique(c_post_B_ni$df_D$time)),
    cue = tar_cue(mode = "thorough")
  )
)     # end target list
