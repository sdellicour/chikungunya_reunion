#!/bin/bash
#SBATCH --job-name=Alignment_101025_glm
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=50000m
#SBATCH --gres="gpu:1"
#SBATCH --time=30-00:00:00

module load beagle-lib/4.0.1-GCC-12.3.0-CUDA-12.1.1

export OMP_NUM_THREADS=8
export MKL_NUM_THREADS=8

java -jar beast_1_10_version_1_10_5_30-12-2024.jar -beagle_gpu -beagle_double -beagle_order 1 -overwrite Alignment_101025_glm.xml
