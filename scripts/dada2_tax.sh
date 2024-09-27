#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 64G
#SBATCH -c 4
#SBATCH -t 12:00:00

Rscript scripts/dada2_tax.r