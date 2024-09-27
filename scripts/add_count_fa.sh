#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 228G
#SBATCH -c 1
#SBATCH -t 4:00:00

Rscript /home/evaes/eDNA/faststorage/Velux/CoastSequence/Autumn/Bac16s/GitHub_repo/scripts/add_count_fa.r