# submit with:
# sbatch slurm/run_mcmc_beta.job

## Load packages 
here::i_am("slurm/run_mcmc_beta.R")
library(here)
library(qs)

library(targets)
library(data.table)
library(BayesianTools)
#library(scales)

#source("_targets_packages.R")
source("R/luc_track.R")
source("R/luct.R")
#set.seed(586)

obs     <- tar_read(c_obs_ni)
pred_ls <- tar_read(c_pred_ls)

## ---- initialise_predictions, eval=recalc, echo=FALSE-------------------------
# Set the time limits over which to calculate and initialise an array for the $B$ matrices.
# for n_t times, there are n_t-1 changes, so n_t should 
# be 1 more than the number of processors requested 
# e.g. for 2015-2019, n_t = 5, n_p = 4
#v_times <- 1750:2019 # 1750:1779
#v_times <- 1950:2019
v_times <- 1950:2020
n_t <- length(v_times)

# read arguments
# which year in v_times to process
i <- as.numeric(commandArgs(trailingOnly = TRUE))[1]
# directory for the output
dir_output <- commandArgs(trailingOnly = TRUE)[2]

#dir_output <- "output"
#i = 69
i_time <- v_times[i+1]
i_time

# this duplicates reordering in combine_blags
# but the additional check is probably a good idea
# as the order determines vector to matrix conversion ordering
obs$dt_B <- dplyr::arrange(obs$dt_B, time, data_source, u_to, u_from)

## ---- estimate_B_by_MCMC_parallel, eval=recalc, echo=FALSE--------------------
# Parallelise over years, finding the posterior distribution of the $B$ matrix by MCMC.
n_iter <- 600000
thin <- 10 # round(max(1, n_iter/20))
# we want three processors to each run one chain, with the 3 internal chains 
n_chains <- 1 # on the same core (DREAMz uses 3 internal chains by default)
n_cores  <- 3 # number of cores to use

# Prior: Half-Normal
# x <- seq(0, 10000, by = 10)
# y <- dnorm(x, mean = 0, sd = 3000)
# plot(x, y, ylim = c(0, max(y)))
prior <- createTruncatedNormalPrior(
   mean = rep(0, n_u^2), 
   sd   = rep(500, n_u^2), 
   lower = rep(0, n_u^2),
   upper = rep(10000, n_u^2))
# # Prior: uniform
# prior <- createUniformPrior(
  # lower = rep(    0, n_u^2), 
  # upper = rep(6000, n_u^2))
  
setUp <- createBayesianSetup(get_loglik, prior = prior, parallel = FALSE)

# get this years observations
dt_B_obs <- obs$dt_B[time == i_time]
dt_G_obs <- obs$dt_G[time == i_time]
dt_L_obs <- obs$dt_L[time == i_time]
dt_A_obs <- obs$dt_A[time == i_time]
dt_D_obs <- obs$dt_D[time == i_time]

# Initial values for chains 
# use LS predictions as a starting point 
res <- try(v_B_ini <- pred_ls$dt_B[time == i_time, area])
# if they have not been calculated, use random values
if (length(res) != n_u^2) v_B_ini <- runif(n_u^2, 0, 100) # random
v_B_ini[is.na(v_B_ini)] <- 0
m_starter <- matrix(rep(v_B_ini, n_cores), nrow = n_cores, byrow = TRUE)
m_starter[2,] <- rep(0, n_u^2)               # all zeroes
m_starter[n_cores,] <- runif(n_u^2, 0, 1000) # random

settings <- list(iterations = n_iter, thin = thin,  
  nrChains = n_chains, startValue = m_starter, 
  burnin  = min(ceiling(n_iter/20), 12000), message = FALSE)
 
# Start cluster with n_cores cores for n chains and export BayesianTools library
cl <- parallel::makeCluster(n_cores)
# set up libraries etc on each worker.
parallel::clusterEvalQ(cl, {
  library(here)
  library(data.table)
  library(BayesianTools)
  source(here("R", "luc_track.R"))
  source(here("R", "luct.R"))
  NULL
})
parallel::clusterExport(cl, c("n_u", "dt_B_obs", "dt_G_obs", 
  "dt_D_obs", "dt_L_obs"))

# for reproducibility
parallel::clusterSetRNGStream(cl, iseed = 586)

# calculate parallel n chains, for each chain the likelihood will be calculated on one core
system.time(
  out <- parallel::parLapply(cl, 1:n_cores, fun = function(X, setUp, settings) 
    runMCMC(setUp, settings, sampler = "DREAMzs"), setUp, settings)
)

Bmap <- MAP(out[[1]])$parametersMAP
#saveRDS(Bmap, file = here("output", paste0("mcmcB_map", i_time, ".rds")))
qsave(Bmap, file = here(dir_output, paste0("mcmcB_map", i_time, ".qs")))
#saveRDS(out,  file = here("output", paste0("mcmcB_", i_time, ".rds")))
qsave(out,  file = here(dir_output, paste0("mcmcB_", i_time, ".qs")))
parallel::stopCluster(cl)

# Combine the chains
out <- createMcmcSamplerList(out)
summary(out) # check for errors

## ----quit_gracefully, eval=recalc, echo=FALSE---------------------------------
# it used to help if scripts exited R without saving when run on LSF job scheduler
# guessing still true with SLURM job scheduler
i_time
quit(save = "no")
