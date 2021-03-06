library(renv)
library(workflowr)
library(targets)
library(tarchetypes)
source("./_targets_packages.R")

# renv
renv::install("sf")
renv::install("corrplot")
renv::install("BayesianTools")
renv::status()
renv::dependencies()
# writes _targets_packages.R
# expurgated version - without clustermq or future, former won't install on POLAR
tar_renv(extras = c("bs4Dash", "gt", "markdown", "pingr", "rstudioapi", "shiny", "shinybusy", "shinyWidgets", "visNetwork"))
renv::status()

# remove.packages("clustermq")
# renv::remove("clustermq")
# clustermq

# targets
tar_manifest(fields = c("format", "memory", "storage", "retrieval"))
tar_manifest(names = starts_with("c_cor"), fields = c("format", "memory", "storage"))
tar_outdated()
# rebuild if needed
system.time(tar_make())

# may need to convert files to unix line endings e.g. 
# find ./analysis/ -name "*.Rmd" -type f -exec dos2unix {} \;
# dos2unix *.R
# dos2unix R/*.R
# dos2unix analysis/*.Rmd

tar_option_get("cue")

tar_glimpse()
# core detail
tar_visnetwork(allow = starts_with("c_"))
# meta detail
tar_visnetwork(allow = starts_with("m_"))
# all detail
tar_visnetwork()
tar_visnetwork(targets_only = TRUE)
tar_visnetwork(allow = starts_with("c_"), targets_only = TRUE)

tar_meta(fields = warnings)
View(tar_meta(fields = path))
tar_path(c_file_CS)

# clean up, in order of completeness
tar_delete()
tar_prune()
tar_destroy()

# some debugging
tar_read(c_file_FC)

obs <- tar_read(c_obs_unc)
df_uncert <-tar_read(c_df_uncert)
tar_load(c_df_A_AgCensus_Sco)
tar_load(c_df_A_AgCensus_Wal)
tar_load(c_df_A_AgCensus_NIr)
tar_load(c_blag_corine)
tar_load(c_blag_iacs)
tar_load(c_blag_lcc)
tar_load(c_blag_lcm)
tar_load(c_obs_unc)
tar_load(c_blag_crome)
tar_load(c(c_blag_AgCensus, c_blag_CS, c_blag_corine, c_blag_fc, c_blag_iacs, c_blag_lcc, c_blag_lcm))


# workflowr
# wflow_build   = render
# wflow_publish = commit .Rmd, render, commit .html
# wflow_git_push = git push
# wflow_git_config(user.name = "peterlevy", user.email = "plevy@ceh.ac.uk")
# wflow_open("analysis/m_CS_plot.Rmd")
wflow_git_config()
wflow_status()
wflow_build() # default is make = TRUE, only where .html doesn't match .Rmd
wflow_build("analysis/m_data_comparison.Rmd")
wflow_build("analysis/index.Rmd")
wflow_build("analysis/m_uqdata.Rmd")
wflow_build("analysis/m_AgCensus_plot.Rmd")
wflow_build("analysis/m_00_status.Rmd") # force the status file to be rendered
# All tracked files that have been modified
wflow_publish(all = TRUE, message = "Adding FC GL targets")
wflow_publish("analysis/m_AgCensus_plot.Rmd", message = "Adding AgCensus plots")
wflow_publish("analysis/m_00_status.Rmd", message = "Adding FC GL targets")
wflow_publish(c("analysis/index.Rmd", "analysis/m_00_status.Rmd"),
  message = "Adding FC GL targets")
wflow_publish("analysis/m_00_status.Rmd", message = "Updating status")
#wflow_view()
wflow_use_github(organization = "NERC-CEH")
wflow_git_push(dry_run = TRUE)
wflow_git_push()


load("../luc_track/data/v_100m.RData", verbose = TRUE)
s_U <- stack("../luc_track/data/st_U_bda_10km.tif")
# s_U <- stack("../luc_track/data/st_U_bda_1km.tif")
# s_U <- stack("../luc_track/data/st_U_bda_1km_Byte_LZW.tif")
system.time(a_B <- getMatrices_fromStack(s_U))
dimnames(a_B) <- list(names_u, names_u, 1990:2019)
a_B <- set_units(a_B, km^2)
str(a_B)
a_B[, , 29]

save(a_B, file = "a_B.rda")
load(file = "a_B.rda", verbose = TRUE)