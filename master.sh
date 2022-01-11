#!/bin/bash
sbatch <<EOT
#!/bin/bash

#SBATCH -J $2
#SBATCH -p general
#SBATCH -o filename_%j.txt
#SBATCH -e filename_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mamalec@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH "--time="$1
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

srun ~/R-4.1.2/build/bin/Rscript --vanilla ~/ccImpute/master.R $2 $3 $4
exit 0
EOT
