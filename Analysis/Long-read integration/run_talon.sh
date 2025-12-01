#!/bin/sh
#SBATCH --partition=multinode
#SBATCH --constraint=[x2680|x2695]
#SBATCH --nodes=2
#SBATCH --ntasks=36
#SBATCH --ntasks-per-core=1
#SBATCH --time=96:00:00
#SBATCH --mem=108g

##sbatch run_talon.sh

module load samtools
module load mamba_install
source myconda
mamba activate TALON

pwd

talon --f config1.csv \
      --db ../myTest32_10s.db \
      --build hg38 \
      --threads 36 \
      --cov 0.95 \
      --identity 0.95 \
      --o output2
