sbatch slurm/run_sample_Upost_en.job
sbatch slurm/run_sample_Upost_sc.job
sbatch slurm/run_sample_Upost_wa.job
sbatch slurm/run_sample_Upost_ni.job
sbatch --account=short4hr slurm/run_sample_Upost_en.job
sbatch --account=short4hr slurm/run_sample_Upost_sc.job
sbatch --account=short4hr slurm/run_sample_Upost_wa.job
sbatch --account=short4hr slurm/run_sample_Upost_ni.job
sbatch --account=short4hr slurm/run_sample_Upost_uk.job
R CMD BATCH --no-restore --no-save "--args  1 10000 en" slurm/sample_Upost.R "output/output_en/console_sample_Upost1.Rout"  > output/output_en/sample_Upost1.out &
R CMD BATCH --no-restore --no-save "--args  1 10000 sc" slurm/sample_Upost.R "output/output_sc/console_sample_Upost1.Rout"  > output/output_sc/sample_Upost1.out &
R CMD BATCH --no-restore --no-save "--args  1 10000 wa" slurm/sample_Upost.R "output/output_wa/console_sample_Upost1.Rout"  > output/output_wa/sample_Upost1.out &
R CMD BATCH --no-restore --no-save "--args  1 10000 ni" slurm/sample_Upost.R "output/output_ni/console_sample_Upost1.Rout"  > output/output_ni/sample_Upost1.out &
R CMD BATCH --no-restore --no-save "--args  1 1000  en" slurm/sample_Upost.R "output/output_en/console_sample_Upost1.Rout"  > output/output_en/sample_Upost1.out &
R CMD BATCH --no-restore --no-save "--args  1 1000  sc" slurm/sample_Upost.R "output/output_sc/console_sample_Upost1.Rout"  > output/output_sc/sample_Upost1.out &
R CMD BATCH --no-restore --no-save "--args  1 1000  wa" slurm/sample_Upost.R "output/output_wa/console_sample_Upost1.Rout"  > output/output_wa/sample_Upost1.out &
R CMD BATCH --no-restore --no-save "--args  1 1000  ni" slurm/sample_Upost.R "output/output_ni/console_sample_Upost1.Rout"  > output/output_ni/sample_Upost1.out &
R CMD BATCH --no-restore --no-save "--args  1 100  en" slurm/sample_Upost.R "output/output_en/console_sample_Upost1.Rout"  > output/output_en/sample_Upost1.out &
R CMD BATCH --no-restore --no-save "--args  2 100  en" slurm/sample_Upost.R "output/output_en/console_sample_Upost2.Rout"  > output/output_en/sample_Upost2.out &
R CMD BATCH --no-restore --no-save "--args  3 100  en" slurm/sample_Upost.R "output/output_en/console_sample_Upost3.Rout"  > output/output_en/sample_Upost3.out &
R CMD BATCH --no-restore --no-save "--args  4 100  en" slurm/sample_Upost.R "output/output_en/console_sample_Upost4.Rout"  > output/output_en/sample_Upost4.out &
R CMD BATCH --no-restore --no-save "--args  1 100  sc" slurm/sample_Upost.R "output/output_sc/console_sample_Upost_100m_1.Rout"  > output/output_sc/sample_Upost_100m_1.out &
R CMD BATCH --no-restore --no-save "--args  2 100  sc" slurm/sample_Upost.R "output/output_sc/console_sample_Upost_100m_2.Rout"  > output/output_sc/sample_Upost_100m_2.out &
R CMD BATCH --no-restore --no-save "--args  3 100  sc" slurm/sample_Upost.R "output/output_sc/console_sample_Upost_100m_3.Rout"  > output/output_sc/sample_Upost_100m_3.out &
R CMD BATCH --no-restore --no-save "--args  4 100  sc" slurm/sample_Upost.R "output/output_sc/console_sample_Upost_100m_4.Rout"  > output/output_sc/sample_Upost_100m_4.out &
R CMD BATCH --no-restore --no-save "--args  1 100  wa" slurm/sample_Upost.R "output/output_wa/console_sample_Upost_100m_1.Rout"  > output/output_wa/sample_Upost_100m_1.out &
R CMD BATCH --no-restore --no-save "--args  2 100  wa" slurm/sample_Upost.R "output/output_wa/console_sample_Upost_100m_2.Rout"  > output/output_wa/sample_Upost_100m_2.out &
R CMD BATCH --no-restore --no-save "--args  3 100  wa" slurm/sample_Upost.R "output/output_wa/console_sample_Upost_100m_3.Rout"  > output/output_wa/sample_Upost_100m_3.out &
R CMD BATCH --no-restore --no-save "--args  4 100  wa" slurm/sample_Upost.R "output/output_wa/console_sample_Upost_100m_4.Rout"  > output/output_wa/sample_Upost_100m_4.out &
R CMD BATCH --no-restore --no-save "--args  1 100  ni" slurm/sample_Upost.R "output/output_ni/console_sample_Upost_100m_1.Rout"  > output/output_ni/sample_Upost_100m_1.out &
R CMD BATCH --no-restore --no-save "--args  2 100  ni" slurm/sample_Upost.R "output/output_ni/console_sample_Upost_100m_2.Rout"  > output/output_ni/sample_Upost_100m_2.out &
R CMD BATCH --no-restore --no-save "--args  3 100  ni" slurm/sample_Upost.R "output/output_ni/console_sample_Upost_100m_3.Rout"  > output/output_ni/sample_Upost_100m_3.out &
R CMD BATCH --no-restore --no-save "--args  4 100  ni" slurm/sample_Upost.R "output/output_ni/console_sample_Upost_100m_4.Rout"  > output/output_ni/sample_Upost_100m_4.out &
