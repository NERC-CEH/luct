R CMD BATCH --no-restore --no-save "--args 13 600000 output/output_en" slurm/run_mcmc_beta_en.R "output/output_en/console13.Rout"  > output/output_en/mcmcb13.out &
R CMD BATCH --no-restore --no-save "--args 26 600000 output/output_en" slurm/run_mcmc_beta_en.R "output/output_en/console26.Rout"  > output/output_en/mcmcb26.out &

R CMD BATCH --no-restore --no-save "--args 63 600000 output/output_sc" slurm/run_mcmc_beta_sc.R "output/output_sc/console63.Rout"  > output/output_sc/mcmcb63.out &
R CMD BATCH --no-restore --no-save "--args 64 600000 output/output_sc" slurm/run_mcmc_beta_sc.R "output/output_sc/console64.Rout"  > output/output_sc/mcmcb64.out &
R CMD BATCH --no-restore --no-save "--args 70 600000 output/output_sc" slurm/run_mcmc_beta_sc.R "output/output_sc/console70.Rout"  > output/output_sc/mcmcb70.out &

R CMD BATCH --no-restore --no-save "--args 104 600000 output/output_wa" slurm/run_mcmc_beta_wa.R "output/output_wa/console104.Rout"  > output/output_wa/mcmcb104.out &

R CMD BATCH --no-restore --no-save "--args  3 output/output_ni" slurm/run_mcmc_beta_ni.R "output/output_ni/console3.Rout"   > output/output_ni/mcmcb3.out &
R CMD BATCH --no-restore --no-save "--args 67 output/output_ni" slurm/run_mcmc_beta_ni.R "output/output_ni/console67.Rout"  > output/output_ni/mcmcb67.out &
