here::i_am("./run_postU.R")
library(renv)
library(here)
library(targets)
library(tarchetypes)
source("./_targets_packages.R")

Sys.setenv(TAR_PROJECT = "postU")
tar_outdated(script = "_targets_postU.R", store = "_targets_postU")
# rebuild if needed
system.time(tar_make(script = "_targets_postU.R", store = "_targets_postU"))

quit(save = "no")