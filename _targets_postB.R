library(targets)
library(tarchetypes) # for tar_knitr_deps_expr()
library(here) # construct file paths relative to project root
library(fs) # file system operations

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
      mode = "thorough", #c("thorough", "always", "never"),
      command   = TRUE,
      depend    = TRUE,
      format    = TRUE,
      iteration = TRUE,
      file      = TRUE)
)

v_times <- 1900:2020
rerun <- TRUE
if (rerun) {
  cue_mode <- "always"
} else {
  cue_mode <- "thorough"
}

c_obs_en <- tar_read(c_obs_en, store = "_targets_wrangler")
c_obs_sc <- tar_read(c_obs_sc, store = "_targets_wrangler")
c_obs_wa <- tar_read(c_obs_wa, store = "_targets_wrangler")
c_obs_ni <- tar_read(c_obs_ni, store = "_targets_wrangler")


# list of target objects
list(

  # Path to MCMC_Beta SLURM job file
  tar_target(
    c_mcmc_job_en,
    fs::path_rel(here("slurm", "run_mcmc_beta_en.job")),
    format = "file"
  ),

  # Path to MCMC_Beta R file
  tar_target(
    c_mcmc_r_en,
    fs::path_rel(here("slurm", "run_mcmc_beta_en.R")),
    format = "file"
  ),

  # Path to MCMC_Beta SLURM job file
  tar_target(
    c_mcmc_job_sc,
    fs::path_rel(here("slurm", "run_mcmc_beta_sc.job")),
    format = "file"
  ),

  # Path to MCMC_Beta SLURM job file
  tar_target(
    c_mcmc_job_wa,
    fs::path_rel(here("slurm", "run_mcmc_beta_wa.job")),
    format = "file"
  ),

  # Path to MCMC_Beta SLURM job file
  tar_target(
    c_mcmc_job_ni,
    fs::path_rel(here("slurm", "run_mcmc_beta_ni.job")),
    format = "file"
  ),

  # Run a SLURM job to estimate Beta by MCMC
  tar_target(
    c_mcmc_out_en,
    run_mcmc_beta_job(c_mcmc_job_en,
      queue_name = "short-serial-4hr",
      dir_output = "output/output_en",
      v_times = v_times, c_obs_en),
    cue = tar_cue(mode = cue_mode)
  ),

  # Run a SLURM job to estimate Beta by MCMC
  tar_target(
    c_mcmc_out_sc,
    run_mcmc_beta_job(c_mcmc_job_sc,
    queue_name = "short-serial-4hr",
    dir_output = "output/output_sc",
      v_times = v_times, c_obs_sc),
    cue = tar_cue(mode = cue_mode)
  ),

  # Run a SLURM job to estimate Beta by MCMC
  tar_target(
    c_mcmc_out_wa,
    run_mcmc_beta_job(c_mcmc_job_wa,
    queue_name = "short-serial-4hr",
    dir_output = "output/output_wa",
      v_times = v_times, c_obs_wa),
    cue = tar_cue(mode = cue_mode)
  ),

  # Run a SLURM job to estimate Beta by MCMC
  tar_target(
    c_mcmc_out_ni,
    run_mcmc_beta_job(c_mcmc_job_ni,
    queue_name = "short-serial-4hr",
    dir_output = "output/output_ni",
      v_times = v_times, c_obs_ni),
    cue = tar_cue(mode = cue_mode)
  )
)     # end target list
