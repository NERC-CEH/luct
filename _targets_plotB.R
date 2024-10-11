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
    "dplyr", "purrr", "units", "data.table", "ggplot2", "stringr", "qs",
    "zoo", "mgcv", "reshape2", "readxl", "tidyr", "sp", "sf", "raster", #"rgeos", "ggforce", "plyr",
    "rgdal", "grid", "spCEH", "scico", "stars", "BayesianTools", "scales"),
  format = "qs",
  cue = tar_cue(
      mode = "always", #c("thorough", "always", "never"),
      command   = TRUE,
      depend    = TRUE,
      format    = TRUE,
      iteration = TRUE,
      file      = TRUE)
)

start_year <- 1900
end_year   <- 2020
fig_start_year <- 1900
fig_end_year   <- 2020

c_obs_en <- tar_read(c_obs_en, store = "_targets_wrangler")
c_obs_sc <- tar_read(c_obs_sc, store = "_targets_wrangler")
c_obs_wa <- tar_read(c_obs_wa, store = "_targets_wrangler")
c_obs_ni <- tar_read(c_obs_ni, store = "_targets_wrangler")

c_blag_lcm_uk <- tar_read(c_blag_lcm_uk, store = "_targets_wrangler")

c_mcmc_out_en <- tar_read(c_mcmc_out_en, store = "_targets_postB")
c_mcmc_out_sc <- tar_read(c_mcmc_out_sc, store = "_targets_postB")
c_mcmc_out_wa <- tar_read(c_mcmc_out_wa, store = "_targets_postB")
c_mcmc_out_ni <- tar_read(c_mcmc_out_ni, store = "_targets_postB")

# list of target objects
list(

  # Plot the results and write summary output
  tar_target(
    c_post_B_en,
    get_post_plots(
      start_year = start_year,
      end_year = end_year,
      dir_output = "output/output_en",
      v_mcmc_fname_Bmap = c_mcmc_out_en,
      fig_start_year = fig_start_year ,
      fig_end_year   = fig_end_year,
      obs_unc = c_obs_en,
      obs_exc = c_obs_en,
      blag_lcm = c_blag_lcm_uk,
      start  = 10000,
      mcmc_diag_plot_year = 2019)
  ),

  # Plot the results and write summary output
  tar_target(
    c_post_B_sc,
    get_post_plots(
      start_year = start_year,
      end_year = end_year,
      dir_output = "output/output_sc",
      v_mcmc_fname_Bmap = c_mcmc_out_sc,
      fig_start_year = fig_start_year ,
      fig_end_year   = fig_end_year,
      obs_unc = c_obs_sc,
      obs_exc = c_obs_sc,
      #v_data_source = c("AgCensus", "MODIS", "CS", "FC", "IACS"),
      blag_lcm = c_blag_lcm_uk,
      start  = 10000,
      mcmc_diag_plot_year = 2019)
  ),

  # Plot the results and write summary output
  tar_target(
    c_post_B_wa,
    get_post_plots(
      start_year = start_year,
      end_year = end_year,
      dir_output = "output/output_wa",
      v_mcmc_fname_Bmap = c_mcmc_out_wa,
      fig_start_year = fig_start_year ,
      fig_end_year   = fig_end_year,
      obs_unc = c_obs_wa,
      obs_exc = c_obs_wa,
      #v_data_source = c("AgCensus", "MODIS", "CS", "FC", "IACS"),
      blag_lcm = c_blag_lcm_uk,
      start  = 10000,
      mcmc_diag_plot_year = 2019)
  ),

  # Plot the results and write summary output
  tar_target(
    c_post_B_ni,
    get_post_plots(
      start_year = start_year,
      end_year = end_year,
      dir_output = "output/output_ni",
      v_mcmc_fname_Bmap = c_mcmc_out_ni,
      fig_start_year = fig_start_year ,
      fig_end_year   = fig_end_year,
      obs_unc = c_obs_ni,
      obs_exc = c_obs_ni,
      #v_data_source = c("AgCensus", "MODIS", "CS", "FC", "IACS"),
      blag_lcm = c_blag_lcm_uk,
      start  = 10000,
      mcmc_diag_plot_year = 2019)
  )
)     # end target list
