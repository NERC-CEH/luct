library(targets)
library(tarchetypes) # for tar_knitr_deps_expr()
library(here) # construct file paths relative to project root
library(fs) # file system operations

# My main computer supports 32 threads, so why not use them?
#
# Uncomment the following code and use tar_make_future(workers = 30L) to make the targets
# in parallel sessions on the local machine.
#
# library(future)
# future::plan("multicore")
# multicore is *not* recommended in Rstudio sessions, but "multisession" doesn't seem to work for me
#
# Unfortunately, in this project you typically only need to recompute a small number of targets
# so the target-level parallelism doesn't help much.
# In this project, it would probably be more helpful to aim for multithreading *within* each target.

# Define custom functions and other global objects in `R/functions.R`.
# All the functions to calculate the targets are defined here.
source(here::here("R", "luc_track.R"))
source(here::here("R", "luct.R"))


# Set target-specific options such as packages.
# It is recommended to load here all the packages required by functions managed only by {targets}.
# Packages required by {workflowr} notebooks should be loaded in those notebooks
# so they can be run manually (rather than exclusively by {targets}).
options(tidyverse.quiet = TRUE)
options(bitmapType='cairo')
tar_option_set(
  packages = c(
    "dplyr", "purrr", "units", "data.table", "ggplot2", "stringr",
    "zoo", "mgcv", "reshape2", "readxl", "tidyr", "sp", "sf", "raster", #"rgeos", "ggforce", "plyr", 
    "rgdal", "grid", "spCEH", "scico", "stars", "BayesianTools"),
  format = "qs",
  cue = tar_cue(
      mode = "thorough", #c("thorough", "always", "never"),
      command   = TRUE,
      depend    = TRUE,
      format    = TRUE,
      iteration = TRUE,
      file      = TRUE)
)

# list of target objects
list(

  # CORE pipeline targets ----

  # Path to raw AgCensus_Eng data file
  tar_target(
    c_file_AgCensus_Eng,
    c(fs::path_rel(here("data-raw/AgCensus/England", 
        "AgCensus_England_ha_1900-2010.csv")),
      fs::path_rel(here("data-raw/AgCensus/England", 
        "AgCensus_England_ha_1983-2019.csv"))),
    format = "file"
  ),

  # Path to raw AgCensus_Sco data file
  tar_target(
    c_file_AgCensus_Sco,
    c(fs::path_rel(here("data-raw/AgCensus/Scotland", 
        "AgCensus_Scotland_ha_1883-2014.csv")),
      fs::path_rel(here("data-raw/AgCensus/Scotland", 
        "AgCensus_Scotland_ha_2009-2019.csv"))),
    format = "file"
  ),

  # Path to raw AgCensus_Wal data file
  tar_target(
    c_file_AgCensus_Wal,
    c(fs::path_rel(here("data-raw/AgCensus/Wales", 
        "AgCensus_Wales_ha_1867-2012.csv")),
      fs::path_rel(here("data-raw/AgCensus/Wales", 
        "AgCensus_Wales_ha_1998-2019.csv"))),
    format = "file"
  ),

  # Path to raw AgCensus_NIr data file
  tar_target(
    c_file_AgCensus_NIr,
    fs::path_rel(here("data-raw/AgCensus/NIreland", 
        "AgCensus_NIreland_ha_1981-2019.csv")),
    format = "file"
  ),

  # Wrangle AgCensus_Eng data
  tar_target(
    c_df_A_AgCensus_Eng,
    wrangle_AgCensus_Eng(c_file_AgCensus_Eng)
  ),

  # Wrangle AgCensus_Sco data
  tar_target(
    c_df_A_AgCensus_Sco,
    wrangle_AgCensus_Sco(c_file_AgCensus_Sco)
  ),

  # Wrangle AgCensus_Wal data
  tar_target(
    c_df_A_AgCensus_Wal,
    wrangle_AgCensus_Wal(c_file_AgCensus_Wal)
  ),

  # Wrangle AgCensus_NIr data
  tar_target(
    c_df_A_AgCensus_NIr,
    wrangle_AgCensus_NIr(c_file_AgCensus_NIr)
  ),

  # Combine AgCensus data
  tar_target(
    c_blag_AgCensus,
    combine_AgCensus(l_df = list(c_df_A_AgCensus_Eng$df_A, 
                                 c_df_A_AgCensus_Sco$df_A, 
                                 c_df_A_AgCensus_Wal$df_A, 
                                 c_df_A_AgCensus_NIr$df_A))
  ),

  # Path to raw CS data file
  tar_target(
    c_file_CS,
    fs::path_rel(here("data-raw/CS", "UK_LUC_matrices_2018i.csv")),
    format = "file"
  ),

  # Wrangle CS data
  tar_target(
    c_blag_CS,
    wrangle_CS(fpath = c_file_CS)
  ),

  # Path to raw FC data fileS
  tar_target(
    c_file_fc,
    c(fs::path_rel(here("data-raw/FC/timeSeries", 
        "forest_planting_byYear_UK.csv")),
      fs::path_rel(here("data-raw/FC/timeSeries", 
        "Deforestation_Areas_for_CEH_1990-2019.xlsx"))),
    format = "file"
  ),

  # Wrangle FC data
  tar_target(
    c_blag_fc,
    wrangle_FC(c_file_fc)
  ),

  # # Path to test SLURM job file
  # tar_target(
    # c_cor_fname_job,
    # fs::path_rel(here("slurm", "runTest.job")),
    # format = "file"
  # ),

  # # Process test data
  # tar_target(
    # c_cor_job,
    # run_corine_job(c_cor_fname_job)
  # ),
      
  # Path to FC R code file
  tar_target(
    c_fc_fname_code,
    fs::path_rel(here("slurm", "process_FC.R")),
    format = "file"
  ),  
  
  # Run a SLURM job to process FC data
  tar_target(
    c_fc_fname_out,
    run_FC_job(c_fc_fname_code),
    format = "file",
    cue = tar_cue(mode = "thorough")
  ),

  # Path to CORINE SLURM job file
  tar_target(
    c_cor_fname_job,
    fs::path_rel(here("slurm", "process_CORINE.job")),
    format = "file"
  ),

  # Run a SLURM job to process CORINE data
  tar_target(
    c_level1_cor,
    run_corine_job(c_cor_fname_job),
    cue = tar_cue(mode = "thorough")
  ),
  
  # Get BLAG from CORINE Level1 output
  tar_target(
    c_blag_corine,
    getBLAG_fromU(
      v_times  = c_level1_cor$v_times,
      v_fnames = c_level1_cor$v_fnames,
      name_data_source = "CORINE", names_u),
    cue = tar_cue(mode = "thorough")
  ),

  # Path to IACS SLURM job file
  tar_target(
    c_iacs_fname_job,
    fs::path_rel(here("slurm", "process_IACS.job")),
    format = "file"
  ),

  # Run a SLURM job to process IACS data
  tar_target(
    c_level1_iacs,
    run_iacs_job(c_iacs_fname_job),
    cue = tar_cue(mode = "thorough")
  ),
    
  # Get BLAG from IACS Level1 output
  tar_target(
    c_blag_iacs,
    getBLAG_fromU(
      v_times  = c_level1_iacs$v_times,
      v_fnames = c_level1_iacs$v_fnames,
      name_data_source = "IACS", names_u),
    cue = tar_cue(mode = "thorough")
  ),
  
  # Path to LCC SLURM job file
  tar_target(
    c_lcc_fname_job,
    fs::path_rel(here("slurm", "process_LCC.job")),
    format = "file"
  ),

  # Run a SLURM job to process LCC data
  tar_target(
    c_level1_lcc,
    run_lcc_job(c_lcc_fname_job),
    cue = tar_cue(mode = "thorough")
  ),
    
  # Get BLAG from LCC Level1 output
  tar_target(
    c_blag_lcc,
    getBLAG_fromU(
      v_times  = c_level1_lcc$v_times,
      v_fnames = c_level1_lcc$v_fnames,
      name_data_source = "LCC", names_u),
    cue = tar_cue(mode = "thorough")
  ),
  
  # Path to LCM SLURM job file
  tar_target(
    c_lcm_fname_job,
    fs::path_rel(here("slurm", "process_LCM.job")),
    format = "file"
  ),

  # Run a SLURM job to process LCM data
  tar_target(
    c_level1_lcm,
    run_lcm_job(c_lcm_fname_job),
    cue = tar_cue(mode = "thorough")
  ),
    
  # Apply a land mask to the LCM data
  tar_target(
    c_level1_lcm_masked,
    apply_mask(c_level1_lcm),
    cue = tar_cue(mode = "thorough")
  ),
    
  # Get BLAG from LCM Level1 output
  tar_target(
    c_blag_lcm,
    getBLAG_fromU(
      v_times  = c_level1_lcm_masked$v_times,
      v_fnames = c_level1_lcm_masked$v_fnames,
      name_data_source = "LCM", names_u),
    cue = tar_cue(mode = "never")
  ),
        
  # Run a SLURM job to process CROME data
  tar_target(
    c_level1_crome,
    run_crome_job(c_crome_fname_job),
    cue = tar_cue(mode = "thorough")
  ),
    
  # Get BLAG from CROME Level1 output
  tar_target(
    c_blag_crome,
    getBLAG_fromU(
      v_times  = c_level1_crome$v_times,
      v_fnames = c_level1_crome$v_fnames,
      name_data_source = "CROME", names_u),
    cue = tar_cue(mode = "thorough")
  ),
        
  # Combine BLAGs to give list of data tables with all obs
  tar_target(
    c_obs_all,
    combine_blags(
      l_blags = list(c_blag_AgCensus, c_blag_CS, c_blag_corine, c_blag_fc, c_blag_iacs, c_blag_lcc, c_blag_lcm, c_blag_crome)),
    cue = tar_cue(mode = "thorough")
  ),
            
  # Exclude some data sources which we do not want to use
  tar_target(
    c_obs,
    set_exclusions(
      c_obs_all),
    cue = tar_cue(mode = "thorough")
  ),
               
  # Calculate the relative uncertainty for the data sources
  tar_target(
    c_df_uncert,
    get_uncert_scaling(
      c_obs, 
      v_names_sources = 
        c("AgCensus", "CS", "FC", "LCM", "CORINE", "LCC", "IACS", "CROME"),
        cv_AgCensus = 0.1)
  ),
                
  # Add the relative uncertainties to the data sources
  tar_target(
    c_obs_unc,
    add_uncert(
      c_obs, 
      c_df_uncert)
  ),
                
  # Predict the Beta matrix by least-squares
  tar_target(
    c_pred_ls,
    get_pred_ls(c_obs, start_year = 1990, end_year = 2019),
    cue = tar_cue(mode = "never")
  ),
                    
  # # Predict the posterior Beta matrix by MCMC in serial
  # tar_target(
    # c_B_post,
    # get_post_mcmc_serial(c_obs, c_pred_ls, start_year = 2017, end_year = 2019, n_iter = 1000)
  # ),
  
  # # Path to MCMC_Beta SLURM job file
  # tar_target(
    # c_mcmc_fname_job,
    # fs::path_rel(here("slurm", "run_mcmc_beta.job")),
    # format = "file"
  # ),

  # # Run a SLURM job to estimate Beta by MCMC
  # tar_target(
    # c_mcmc_fname_out,
    # run_mcmc_beta(c_mcmc_fname_job),
    # cue = tar_cue(mode = "thorough")
  # ),    
  
  # # META pipeline targets ----

  # # Any META steps that are straight computation (not Rmd notebooks)
  # # are treated exactly like CORE targets.
  # #
  # # Each META pipeline should have leaves that are the rendering
  # # of one or more {workflowr} Rmd notebooks.
  # # These notebooks are part-managed by {targets}
  # # and part-managed by {workflowr}.
  # #
  # # For the moment, I only get {targets} to monitor the input Rmd files
  # # and any data dependencies (via tar_read and tar_load).
  # # That is, I do NOT include the rendered HTML outputs
  # # as files to be watched by {targets}.
  # # This means I am only relying on {targets} to rebuild the reports if
  # # any upstream dependencies, including the Rmd file, change.
  # #
  # # The targets defined here use workflowr::wflow_build()
  # # to render the Rmd file.
  # # Because the target is not monitoring the rendered HTML output
  # # it can't tell whether the report has been rendered.
  # # However, {workflowr} *does* monitor the relationship
  # # between the Rmd and HTML file,
  # # so wflow_status() knows whether the rendering is up to date.
  # # (If I included the rendered HTML report as a dependency
  # # targets would flag the report as invalidated if I manually
  # # rebuilt or published the report with workflowr.)
  # #
  # # The user needs to manually execute wflow_publish()
  # # to publish the report.
  # # NOTE when {workflowr} is run manually (not via {targets})
  # # the usual setup done by _targets.R (e.g. loading packages)
  # # *won't* have been done.
  # # So the {workflowr} Rmd notebooks must do all their own setup
  # # on the assumption they are run manually.
  # #
  # # I should really develop a tarchetype tar_wflow_build like tar_render,
  # # but I don't feel up to that yet,
  # # so I will take the more long-winded tar_target approach

  # # The special {workflowr} Rmd files (index.Rmd, about.Rmd, license.Rmd)
  # # don't generally need to have targets defined here
  # # because they don't generally have any upstream dependencies.
  # # I am relying entirely on {workflowr} to track these files.

  ## m_01 Plot the AgCensus data ----

  # Report plotting the AgCensus data
  tar_target(
    m_AgCensus_plot,
    command = {
      # Scan for targets of tar_read() and tar_load()
      !!tar_knitr_deps_expr(here("analysis", "m_AgCensus_plot.Rmd"))

      # Build the report
      workflowr::wflow_build(
        here("analysis", "m_AgCensus_plot.Rmd")
      )

      # Track the input Rmd file (and not the rendered HTML file).
      # Make the path relative to keep the project portable.
      fs::path_rel(here("analysis", "m_AgCensus_plot.Rmd"))
    },
    # Track the files returned by the command
    format = "file",
    cue = tar_cue(mode = "thorough")
  ),   # end m_AgCensus_plot
  
  ## m_02 Plot the CS data ----

  # Report investigating how to read the raw CS data
  tar_target(
    m_CS_plot,
    command = {
      # Scan for targets of tar_read() and tar_load()
      !!tar_knitr_deps_expr(here("analysis", "m_CS_plot.Rmd"))

      # Build the report
      workflowr::wflow_build(
        here("analysis", "m_CS_plot.Rmd")
      )

      # Track the input Rmd file (and not the rendered HTML file).
      # Make the path relative to keep the project portable.
      fs::path_rel(here("analysis", "m_CS_plot.Rmd"))
    },
    # Track the files returned by the command
    format = "file",
    cue = tar_cue(mode = "thorough")
  ),   # end m_CS_plot
  
  ## m_03 Plot the data comparison ----

  # Report investigating how to read the raw CS data
  tar_target(
    m_data_comparison,
    command = {
      # Scan for targets of tar_read() and tar_load()
      !!tar_knitr_deps_expr(here("analysis", "m_data_comparison.Rmd"))

      # Build the report
      workflowr::wflow_build(
        here("analysis", "m_data_comparison.Rmd")
      )

      # Track the input Rmd file (and not the rendered HTML file).
      # Make the path relative to keep the project portable.
      fs::path_rel(here("analysis", "m_data_comparison.Rmd"))
    },
    # Track the files returned by the command
    format = "file",
    cue = tar_cue(mode = "thorough")
  ),   # end m_data_comparison
  
  ## m_04 Quantify relative uncertainties ----

  # Report investigating how to read the raw CS data
  tar_target(
    m_uqdata,
    command = {
      # Scan for targets of tar_read() and tar_load()
      !!tar_knitr_deps_expr(here("analysis", "m_uqdata.Rmd"))

      # Build the report
      workflowr::wflow_build(
        here("analysis", "m_uqdata.Rmd")
      )

      # Track the input Rmd file (and not the rendered HTML file).
      # Make the path relative to keep the project portable.
      fs::path_rel(here("analysis", "m_uqdata.Rmd"))
    },
    # Track the files returned by the command
    format = "file",
    cue = tar_cue(mode = "thorough")
  )   # end m_uqdata
)     # end target list
