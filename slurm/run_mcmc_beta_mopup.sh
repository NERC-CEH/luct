R CMD BATCH --no-restore --no-save "--args 27 output/output_en" slurm/run_mcmc_beta_en.R "output/output_en/console27.Rout"  > output/output_en/mcmcb27.out &
R CMD BATCH --no-restore --no-save "--args 67 output/output_en" slurm/run_mcmc_beta_en.R "output/output_en/console67.Rout"  > output/output_en/mcmcb67.out &

R CMD BATCH --no-restore --no-save "--args 11 output/output_sc" slurm/run_mcmc_beta_sc.R "output/output_sc/console11.Rout"  > output/output_sc/mcmcb11.out &
R CMD BATCH --no-restore --no-save "--args 35 output/output_sc" slurm/run_mcmc_beta_sc.R "output/output_sc/console35.Rout"  > output/output_sc/mcmcb35.out &

R CMD BATCH --no-restore --no-save "--args 14 output/output_wa" slurm/run_mcmc_beta_wa.R "output/output_wa/console14.Rout"  > output/output_wa/mcmcb14.out &

R CMD BATCH --no-restore --no-save "--args  3 output/output_ni" slurm/run_mcmc_beta_ni.R "output/output_ni/console3.Rout"   > output/output_ni/mcmcb3.out &
R CMD BATCH --no-restore --no-save "--args 67 output/output_ni" slurm/run_mcmc_beta_ni.R "output/output_ni/console67.Rout"  > output/output_ni/mcmcb67.out &
