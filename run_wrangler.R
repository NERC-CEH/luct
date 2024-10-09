here::i_am("./run_wrangler.R")
library(renv)
library(here)
library(targets)
library(tarchetypes)
source("./_targets_packages.R")

Sys.setenv(TAR_PROJECT = "wrangler")
tar_outdated(script = "_targets_wrangler.R", store = "_targets_wrangler")
# rebuild if needed
system.time(tar_make(script = "_targets_wrangler.R", store = "_targets_wrangler"))

tar_load(c_file_AgCensus_Eng, store = "_targets_wrangler")
quit(save = "no")
