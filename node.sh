#!/bin/bash 
#SBATCH -J opx
#SBATCH -o 2dising%j.txt
##SBATCH --nodelist=compute-0-5
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --mem 230000
##SBATCH -t 13-00:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=schang21@Central.uh.edu
module load intel
module load gcc/6.4.0
module load python/3.7.3

python oopex.py kryovdiction.py kryovmultiply.py kryovlinear_algebra.py kryovobservable.py
