library(targets)
library(tarchetypes) # for tar_knitr_deps_expr()
library(here) # construct file paths relative to project root
library(fs) # file system operations

source(here::here("R", "luc_track.R"))
source(here::here("R", "luct.R"))
source(here::here("R", "luct_AgCensus.R"))

# Set target-specific options such as packages.
# It is recommended to load here all the packages required by functions managed only by {targets}.
# Packages required by {workflowr} notebooks should be loaded in those notebooks
# so they can be run manually (rather than exclusively by {targets}).
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

start_year <- 1900
end_year <- 2020

# list of target objects
list(

  # CORE pipeline targets ----

  # Path to historical UK AgCensus data file
  tar_target(
    c_file_agc_hist,
    fs::path_rel(here("data-raw/AgCensus",
      "all_LU_stats_FULLSERIES_km2.csv")),
    format = "file"
  ),

  # Path to raw AgCensus_Eng data file
  tar_target(
    c_file_AgCensus_Eng,
      fs::path_rel(here("data-raw/AgCensus/England",
        "AgCensus_England_ha_1983-2020.csv")),
    format = "file"
  ),

  # Path to raw AgCensus_Sco data file
  tar_target(
    c_file_AgCensus_Sco,
    c(fs::path_rel(here("data-raw/AgCensus/Scotland",
        "AgCensus_Scotland_ha_1883-2014.csv")),
      fs::path_rel(here("data-raw/AgCensus/Scotland",
        "AgCensus_Scotland_ha_2009-2020.csv"))),
    format = "file"
  ),

  # Path to raw AgCensus_Wal data file
  tar_target(
    c_file_AgCensus_Wal,
    c(fs::path_rel(here("data-raw/AgCensus/Wales",
        "AgCensus_Wales_ha_1867-2012.csv")),
      fs::path_rel(here("data-raw/AgCensus/Wales",
        "AgCensus_Wales_ha_1998-2020.csv"))),
    format = "file"
  ),

  # Path to raw AgCensus_NIr data file
  tar_target(
    c_file_AgCensus_NIr,
    fs::path_rel(here("data-raw/AgCensus/NIreland",
        "AgCensus_NIreland_ha_1981-2020.csv")),
    format = "file"
  ),

  # Wrangle historical AgCensus data
  tar_target(
    c_agc_hist,
    wrangle_AgCensus_historical(c_file_agc_hist)
  ),

  # Process AgCensus_Eng data
  tar_target(
    c_agc_en,
    process_AgCensus_Eng(c_file_AgCensus_Eng, c_agc_hist)
  ),

  # Process AgCensus_Sco data
  tar_target(
    c_agc_sc,
    process_AgCensus_Sco(c_file_AgCensus_Sco, c_agc_hist)
  ),

  # Process AgCensus_Wal data
  tar_target(
    c_agc_wa,
    process_AgCensus_Wal(c_file_AgCensus_Wal, c_agc_hist)
  ),

  # Process AgCensus_NIr data
  tar_target(
    c_agc_ni,
    process_AgCensus_NIr(c_file_AgCensus_NIr, c_agc_hist)
  ),

  # Combine AgCensus data
  tar_target(
    c_blag_AgCensus,
    combine_AgCensus(l_df = list(c_agc_en$df_long,
                                 c_agc_sc$df_long,
                                 c_agc_wa$df_long,
                                 c_agc_ni$df_long),
                                 v_region = v_region)
  ),

  # # Combine AgCensus data using historical data only
  # tar_target(
  #   c_blag_AgCensus,
  #   combine_AgCensus(l_df = list(c_agc_hist$df_long),
  #                                v_region = v_region)
  # ),

  # Path to MODIS data file
  tar_target(
    c_file_MODIS,
    fs::path_rel(here("data-raw/MODIS", "FAOStat_Artificial_land_surface_MODIS.csv")),
    format = "file"
  ),

  # Wrangle MODIS data
  tar_target(
    c_blag_MODIS,
    wrangle_MODIS_urban(fpath = c_file_MODIS)
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
        "forest_planting_byYear_ha.csv")),
      fs::path_rel(here("data-raw/FC/timeSeries",
        "Deforestation_Areas_for_CEH.xlsx"))),
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
    c_blag_corine_uk,
    getBLAG_fromU(
      v_times  = c_level1_cor$v_times,
      v_fnames = c_level1_cor$v_fnames,
      name_data_source = "CORINE",
      region = "uk", names_u),
    cue = tar_cue(mode = "thorough")
  ),

  # Get BLAG from CORINE Level1 output
  tar_target(
    c_blag_corine_en,
    getBLAG_fromU(
      v_times  = c_level1_cor$v_times,
      v_fnames = c_level1_cor$v_fnames,
      name_data_source = "CORINE",
      region = "en", names_u),
    cue = tar_cue(mode = "thorough")
  ),

  # Get BLAG from CORINE Level1 output
  tar_target(
    c_blag_corine_sc,
    getBLAG_fromU(
      v_times  = c_level1_cor$v_times,
      v_fnames = c_level1_cor$v_fnames,
      name_data_source = "CORINE",
      region = "sc", names_u),
    cue = tar_cue(mode = "thorough")
  ),

  # Get BLAG from CORINE Level1 output
  tar_target(
    c_blag_corine_wa,
    getBLAG_fromU(
      v_times  = c_level1_cor$v_times,
      v_fnames = c_level1_cor$v_fnames,
      name_data_source = "CORINE",
      region = "wa", names_u),
    cue = tar_cue(mode = "thorough")
  ),

  # Get BLAG from CORINE Level1 output
  tar_target(
    c_blag_corine_ni,
    getBLAG_fromU(
      v_times  = c_level1_cor$v_times,
      v_fnames = c_level1_cor$v_fnames,
      name_data_source = "CORINE",
      region = "ni", names_u),
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
    c_blag_iacs_en,
    getBLAG_fromU(
      v_times  = c_level1_iacs$v_times,
      v_fnames = c_level1_iacs$v_fnames,
      name_data_source = "IACS",
      region = "en", names_u),
    cue = tar_cue(mode = "never")
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
    c_blag_lcc_en,
    getBLAG_fromU(
      v_times  = c_level1_lcc$v_times,
      v_fnames = c_level1_lcc$v_fnames,
      name_data_source = "LCC",
      region = "en", names_u),
    cue = tar_cue(mode = "never")
  ),

  # Get BLAG from LCC Level1 output
  tar_target(
    c_blag_lcc_sc,
    getBLAG_fromU(
      v_times  = c_level1_lcc$v_times,
      v_fnames = c_level1_lcc$v_fnames,
      name_data_source = "LCC",
      region = "sc", names_u),
    cue = tar_cue(mode = "never")
  ),

  # Get BLAG from LCC Level1 output
  tar_target(
    c_blag_lcc_wa,
    getBLAG_fromU(
      v_times  = c_level1_lcc$v_times,
      v_fnames = c_level1_lcc$v_fnames,
      name_data_source = "LCC",
      region = "wa", names_u),
    cue = tar_cue(mode = "never")
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
    cue = tar_cue(mode = "never")
  ),

  # Get BLAG from LCM Level1 output
  tar_target(
    c_blag_lcm_en,
    getBLAG_fromU(
      v_times  = c_level1_lcm_masked$v_times,
      v_fnames = c_level1_lcm_masked$v_fnames,
      name_data_source = "LCM",
      region = "en", names_u),
    cue = tar_cue(mode = "thorough")
  ),

  # Get BLAG from LCM Level1 output
  tar_target(
    c_blag_lcm_sc,
    getBLAG_fromU(
      v_times  = c_level1_lcm_masked$v_times,
      v_fnames = c_level1_lcm_masked$v_fnames,
      name_data_source = "LCM",
      region = "sc", names_u),
    cue = tar_cue(mode = "thorough")
  ),

  # Get BLAG from LCM Level1 output
  tar_target(
    c_blag_lcm_wa,
    getBLAG_fromU(
      v_times  = c_level1_lcm_masked$v_times,
      v_fnames = c_level1_lcm_masked$v_fnames,
      name_data_source = "LCM",
      region = "wa", names_u),
    cue = tar_cue(mode = "thorough")
  ),

  # Get BLAG from LCM Level1 output
  tar_target(
    c_blag_lcm_ni,
    getBLAG_fromU(
      v_times  = c_level1_lcm_masked$v_times,
      v_fnames = c_level1_lcm_masked$v_fnames,
      name_data_source = "LCM",
      region = "ni", names_u),
    cue = tar_cue(mode = "thorough")
  ),

  # Get BLAG from LCM Level1 output
  tar_target(
    c_blag_lcm_uk,
    getBLAG_fromU(
      v_times  = c_level1_lcm_masked$v_times,
      v_fnames = c_level1_lcm_masked$v_fnames,
      name_data_source = "LCM",
      region = "uk", names_u),
    cue = tar_cue(mode = "thorough")
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
      name_data_source = "CROME",
      region = "en", names_u),
    cue = tar_cue(mode = "never")
  ),

  # Combine BLAGs to give list of data tables with all obs
  tar_target(
    c_obs_all,
    combine_blags(
      l_blags = list(c_blag_AgCensus, c_blag_MODIS, c_blag_CS, c_blag_fc,
      c_blag_corine_uk,
      c_blag_corine_en,
      c_blag_corine_sc,
      c_blag_corine_wa,
      c_blag_corine_ni,
      c_blag_crome, c_blag_iacs_en,
      c_blag_lcc_en,
      c_blag_lcc_sc,
      c_blag_lcc_wa,
      c_blag_lcm_en,
      c_blag_lcm_sc,
      c_blag_lcm_wa,
      c_blag_lcm_ni,
      c_blag_lcm_uk)),
    cue = tar_cue(mode = "thorough")
  ),

  # Convert area in BLAGs from m2 to km2
  tar_target(
    c_obs_km2,
    convert_units(
      c_obs_all,
      old_unit = "m^2",
      new_unit = "km^2"),
    cue = tar_cue(mode = "thorough")
  ),

  # remove this target to meta notebook - one-off operation with all data sets
  # we now want flexibility to exclude some data sources, based on the results from all
  # # Calculate the relative uncertainty for the data sources
  # tar_target(
    # c_df_uncert,
    # get_uncert_scaling(
      # c_obs_exc,
        # v_names_sources =
        # c("AgCensus", "MODIS", "CS", "FC", "LCM", "CORINE", "LCC", "IACS", "CROME"),
        # v_interval_length_sources =
        # c( 1,          1,       8,    1,    3,     6,        1,     1,      1),
        # # this gives a way of removing effect of data sources: long int_lth means effects is reduced by sqrt(n)
        # #c( 1,          1,       8,    1,    1e9,   1e9,      1e9,   1,    1e9),
        # v_start_year_source =
        # c( 1750,       2001,    1950, 1900, 1990,  2000,     2015,  2004,   2016),
        # cv_AgCensus = 0.1),
    # cue = tar_cue(mode = "thorough")
  # ),

  # replacing above, this target just holds a data file of previously-calculated uncertainties
  # we now want flexibility to exclude some data sources, based on the results from all
  # Read in the relative uncertainty for the data sources
  tar_target(
    c_fname_df_uncert,
      fs::path_rel(here("data", "df_uncert.qs")),
    format = "file"
  ),

  # Interpolate NAs in BLAGs
  tar_target(
    c_obs_filled,
    interpolate_blag(
      c_obs_km2,
      start_year = start_year,
      end_year   = end_year),
    cue = tar_cue(mode = "thorough")
  ),

  # Add the relative uncertainties to the data sources
  tar_target(
    c_obs_unc,
    add_uncert(
      c_obs_filled,
      c_fname_df_uncert),
    cue = tar_cue(mode = "thorough")
  ),

  # Exclude some data sources which we do not want to use
  tar_target(
    c_obs_en,
    set_exclusions(
      c_obs_unc,
      regions_toInclude = "en"),
    cue = tar_cue(mode = "thorough")
  ),

  # Exclude some data sources which we do not want to use
  tar_target(
    c_obs_sc,
    set_exclusions(
      c_obs_unc,
      regions_toInclude = "sc"),
    cue = tar_cue(mode = "thorough")
  ),

  # Exclude some data sources which we do not want to use
  tar_target(
    c_obs_wa,
    set_exclusions(
      c_obs_unc,
      regions_toInclude = "wa"),
    cue = tar_cue(mode = "thorough")
  ),

  # Exclude some data sources which we do not want to use
  tar_target(
    c_obs_ni,
    set_exclusions(
      c_obs_unc,
      regions_toInclude = "ni"),
    cue = tar_cue(mode = "thorough")
  ),

  # Plot the observations
  tar_target(
    c_obs_plots_en,
    get_obs_plots(
      start_year = start_year,
      end_year = end_year,
      dir_output = "output/output_en",
      fig_start_year = start_year ,
      fig_end_year   = end_year,
      obs_unc = c_obs_en,
      obs_exc = c_obs_en)
  ),

  # Plot the observations
  tar_target(
    c_obs_plots_sc,
    get_obs_plots(
      start_year = start_year,
      end_year = end_year,
      dir_output = "output/output_sc",
      fig_start_year = start_year ,
      fig_end_year   = end_year,
      obs_unc = c_obs_sc,
      obs_exc = c_obs_sc)
  ),

  # Plot the observations
  tar_target(
    c_obs_plots_wa,
    get_obs_plots(
      start_year = start_year,
      end_year = end_year,
      dir_output = "output/output_wa",
      fig_start_year = start_year ,
      fig_end_year   = end_year,
      obs_unc = c_obs_wa,
      obs_exc = c_obs_wa)
  ),

  # Plot the observations
  tar_target(
    c_obs_plots_ni,
    get_obs_plots(
      start_year = start_year,
      end_year = end_year,
      dir_output = "output/output_ni",
      fig_start_year = start_year ,
      fig_end_year   = end_year,
      obs_unc = c_obs_ni,
      obs_exc = c_obs_ni)
  )
)     # end target list
