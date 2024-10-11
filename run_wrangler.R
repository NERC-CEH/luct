here::i_am("./run_wrangler.R")
library(renv)
library(here)
library(targets)
library(tarchetypes)
source("./_targets_packages.R")

Sys.setenv(TAR_PROJECT = "wrangler")
tar_outdated(script = "_targets_wrangler.R", store = "_targets_wrangler")
# run in serial
system.time(tar_make(script = "_targets_wrangler.R", store = "_targets_wrangler"))
# run in parallel
system.time(tar_make_future(workers = 4L, script = "_targets_wrangler.R", store = "_targets_wrangler"))

obs_exc <- tar_read(c_obs_en, store = "_targets_wrangler")
unique(obs_exc$dt_D$data_source)
dt <- qread(c_fname_df_uncert)
fwrite(dt, file = "data/df_uncert.csv")

quit(save = "no")
