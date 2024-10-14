library(targets)
library(tarchetypes) # for tar_knitr_deps_expr()
library(here) # construct file paths relative to project root
library(fs) # file system operations

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

res <- 10000
n_samples <- 4

rerun <- TRUE
if (rerun) {
  cue_mode <- "always"
} else {
  cue_mode <- "thorough"
}

# list of target objects
list(

  # Write a slurm job file for England
  tar_target(
    fname_job_en,
    write_postU_job(path_name = here("slurm"),
      region = "en",
      res = res,
      n_samples = n_samples)
  ),
  # Write a slurm job file for Scotland
  tar_target(
    fname_job_sc,
    write_postU_job(path_name = here("slurm"),
      region = "sc",
      res = res,
      n_samples = n_samples)
  ),
  # Write a slurm job file for Wales
  tar_target(
    fname_job_wa,
    write_postU_job(path_name = here("slurm"),
      region = "wa",
      res = res,
      n_samples = n_samples)
  ),
  # Write a slurm job file for N Ireland
  tar_target(
    fname_job_ni,
    write_postU_job(path_name = here("slurm"),
      region = "ni",
      res = res,
      n_samples = n_samples)
  ),

  # Run a SLURM job to sample posterior U from Beta
  tar_target(
    c_fname_dt_Umap_en,
    command = {
      # Scan for targets of tar_read() and tar_load()
      !!tar_knitr_deps_expr(here("slurm", "postU.R"))
      fs::path_rel(here("slurm", "postU.R"))
      sample_Upost_job(fname_job_en, region = "en", res = res)
    },
    # Track the files returned by the command
    format = "file",
    cue = tar_cue(mode = cue_mode)
  ),

  # Run a SLURM job to sample posterior U from Beta
  tar_target(
    c_fname_dt_Umap_sc,
    command = {
      # Scan for targets of tar_read() and tar_load()
      !!tar_knitr_deps_expr(here("slurm", "postU.R"))
      fs::path_rel(here("slurm", "postU.R"))
      sample_Upost_job(fname_job_sc, region = "sc", res = res)
    },
    # Track the files returned by the command
    format = "file",
    cue = tar_cue(mode = cue_mode)
  ),

  # Run a SLURM job to sample posterior U from Beta
  tar_target(
    c_fname_dt_Umap_wa,
    command = {
      # Scan for targets of tar_read() and tar_load()
      !!tar_knitr_deps_expr(here("slurm", "postU.R"))
      fs::path_rel(here("slurm", "postU.R"))
      sample_Upost_job(fname_job_wa, region = "wa", res = res)
    },
    # Track the files returned by the command
    format = "file",
    cue = tar_cue(mode = cue_mode)
  ),

  # Run a SLURM job to sample posterior U from Beta
  tar_target(
    c_fname_dt_Umap_ni,
    command = {
      # Scan for targets of tar_read() and tar_load()
      !!tar_knitr_deps_expr(here("slurm", "postU.R"))
      fs::path_rel(here("slurm", "postU.R"))
      sample_Upost_job(fname_job_ni, region = "ni", res = res)
    },
    # Track the files returned by the command
    format = "file",
    cue = tar_cue(mode = cue_mode)
  )
)     # end target list
