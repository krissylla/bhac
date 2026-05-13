#!/bin/bash
#SBATCH --job-name="FrozenSlab"
#SBATCH --nodes=16
#SBATCH --ntasks=16
#SBATCH --time=00:10:00
#SBATCH --partition=rome

module purge
module load 2023
module load OpenMPI/4.1.5-GCC-12.3.0


cd $HOME/bhac_runs/OnlyFrozenSlab

echo "setting up grid"

$BHAC_DIR/setup.pl -d=33 -phi=3 -z=2 -g=12,12,12 -p=rmhd -eos=gamma -nf=0 -arch=gfortran10 -coord=cart

echo "grid set up complete. First make clean:"

make clean
echo "now  we compile"

make



