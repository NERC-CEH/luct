library(renv)
library(workflowr)
library(targets)
library(tarchetypes)

# renv
renv::install("sf")
renv::install("fasterize")
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
# rebuild if needed
system.time(tar_make())
tar_meta(fields = warnings)
View(tar_meta(fields = path))
tar_path(c_file_CS)

# clean up, in order of completeness
tar_delete()
tar_prune()
tar_destroy()

# some debugging
tar_read(c_file_FC)

tar_load(c_df_A_AgCensus_Eng)
tar_load(c_df_A_AgCensus_Sco)
tar_load(c_df_A_AgCensus_Wal)
tar_load(c_df_A_AgCensus_NIr)
tar_load(c_blag_corine)
tar_load(c_blag_iacs)
tar_load(c_blag_lcc)
tar_load(c_blag_lcm)
tar_load(c_blag_CS)
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
wflow_build("analysis/index.Rmd")
wflow_build("analysis/m_uqdata.Rmd")
wflow_build("analysis/m_CS_plot.Rmd")
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
