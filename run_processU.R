here::i_am("./run_processU.R")
library(renv)
library(here)
library(targets)
library(tarchetypes)
source("./_targets_packages.R")

Sys.setenv(TAR_PROJECT = "processU")
tar_outdated(script = "_targets_processU.R", store = "_targets_processU")
# run in serial
system.time(tar_make(script = "_targets_processU.R", store = "_targets_processU"))
# run in parallel
system.time(tar_make_future(workers = 4L, script = "_targets_processU.R", store = "_targets_processU"))

quit(save = "no")
