here::i_am("./run_postB.R")
library(renv)
library(here)
library(targets)
library(tarchetypes)
source("./_targets_packages.R")

Sys.setenv(TAR_PROJECT = "postB")
tar_outdated(script = "_targets_postB.R", store = "_targets_postB")
# rebuild if needed
system.time(tar_make(script = "_targets_postB.R", store = "_targets_postB"))

quit(save = "no")
