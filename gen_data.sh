#!/bin/bash

#SBATCH -J gen_data
#SBATCH -p general
#SBATCH -o filename_%j.txt
#SBATCH -e filename_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mamalec@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH --mem=128G
#SBATCH --qos=highprio

export OMP_NUM_THREADS=16
PATH="$PATH:$HOME/.local/bin"
module unload gcc
module load intel/19.0.5
module load gcc/9.3.0
module load boost/gnu/1.72.0
module load hdf5/intel/serial/1.12.0
module load python

srun ~/R-4.1.2/build/bin/Rscript --vanilla ~/ccImpute/gen_synth_data.R

