#!/bin/bash
#SBATCH --job-name="JetPropagation"
#SBATCH --nodes=16
#SBATCH --ntasks=16
#SBATCH --time=00:600:00
#SBATCH --partition=rome

module purge
module load 2023
module load OpenMPI/4.1.5-GCC-12.3.0




cd $HOME/bhac_runs/JetPropagation
echo "Running the simulation"

mpiexec  ./bhac -i amrvac.par >output/out


