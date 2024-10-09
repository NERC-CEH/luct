here::i_am("./run_plotB.R")
library(renv)
library(here)
library(targets)
library(tarchetypes)
source("./_targets_packages.R")

Sys.setenv(TAR_PROJECT = "plotB")
tar_outdated(script = "_targets_plotB.R", store = "_targets_plotB")
# run in serial
system.time(tar_make(script = "_targets_plotB.R", store = "_targets_plotB"))
# run in parallel
system.time(tar_make_future(workers = 4L, script = "_targets_plotB.R", store = "_targets_plotB"))

quit(save = "no")
