#!/bin/bash -l
#SBATCH -p batch
#SBATCH -c 8
#SBATCH --time=0-00:30:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tim.seure@uni.lu

module load compiler/Go/1.22.1
go run main.go