#!/bin/sh

#SBATCH --job-name=rb_est
#SBATCH --account=pi-xinhe
#SBATCH --partition=caslake
#SBATCH --ntasks-per-node=1 # number of tasks
#SBATCH --cpus-per-task=2 # number of threads
#SBATCH --mem=168GB
#SBATCH --time=36:00:00
#SBATCH --output=/project/xinhe/lifan/neuron_stim/rb_est_sbatch.log

module load gsl/2.6
module load R/4.2.0

Rscript /project/xinhe/lifan/neuron_stim/mateqtl_100lines_output/rb_est/rb_ns_vs_gtex.R