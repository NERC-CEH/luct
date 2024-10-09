library(targets)
library(tarchetypes) # for tar_knitr_deps_expr()
library(here) # construct file paths relative to project root
library(fs) # file system operations
library(future)
future::plan("multicore")

source(here::here("R", "luc_track.R"))
source(here::here("R", "luct.R"))
source(here::here("R", "luct_AgCensus.R"))

# Set target-specific options such as packages.
# It is recommended to load here all the packages required by functions managed only by {targets}.
# Packages required by {workflowr} notebooks should be loaded in those notebooks
# so they can be run manually (rather than exclusively by {targets}).
options(tidyverse.quiet = TRUE)
options(bitmapType='cairo')
tar_option_set(
  packages = c(
    "dplyr", "purrr", "units", "data.table", "ggplot2", "stringr", "qs", "ggthemes",
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

v_region <- c("en", "sc", "wa", "ni")
v_dir_output <- paste0("output/output_", v_region)
res <- 10000
i_sample <- 1

# list of target objects
list(

# Read postU output
  tar_target(
    v_fname_dt_luv,
    paste0(v_dir_output, "/dt_luv_", v_region, "_smp", i_sample, "_", res, "m.qs"),
    format = "file"
  ),

  tar_target(
    v_fname_dt_lu,
    paste0(v_dir_output, "/dt_u_", v_region, "_smp", i_sample, "_", res, "m.qs"),
    format = "file"
  ),

  tar_target(
    l_dt_luv,
    lapply(v_fname_dt_luv, qread)
  ),

  tar_target(
    l_dt_lu,
    lapply(v_fname_dt_lu, qread)
  ),

  tar_target(
    l_p_plotU_en,
    plot_postU(
      l_dt_luv[[1]],
      l_dt_lu[[1]],
      n_vectors_to_plot = 1000,
      times_achangin = c(1901, 2019)
    )
  ),
  tar_target(
    l_p_plotU_sc,
    plot_postU(
      l_dt_luv[[2]],
      l_dt_lu[[2]],
      n_vectors_to_plot = 1000,
      times_achangin = c(1901, 2019)
    )
  ),
  tar_target(
    l_p_plotU_wa,
    plot_postU(
      l_dt_luv[[3]],
      l_dt_lu[[3]],
      n_vectors_to_plot = 1000,
      times_achangin = c(1901, 2019)
    )
  ),
  tar_target(
    l_p_plotU_ni,
    plot_postU(
      l_dt_luv[[4]],
      l_dt_lu[[4]],
      n_vectors_to_plot = 1000,
      times_achangin = c(1901, 2019)
    )
  )
)     # end target list
