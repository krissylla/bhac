#!/bin/bash
#SBATCH --job-name="FrozenSlab"
#SBATCH --nodes=16
#SBATCH --ntasks=16
#SBATCH --time=00:300:00
#SBATCH --partition=rome

module purge
module load 2023
module load OpenMPI/4.1.5-GCC-12.3.0

echo 'export BHAC_DIR=$HOME/codes/bhac' >> ~/.bashrc
source ~/.bashrc

cp –r ~/bhac_runs/FrozenSlab $TMPDIR/bhac_runs/.
cp -r $BHAC_DIR $TMPDIR/codes/bhac
cd $TMPDIR/bhac_runs/FrozenSlab

echo 'export BHAC_DIR=$HOME/codes/bhac' >> ~/.bashrc
source ~/.bashrc

echo "setting up grid"

$BHAC_DIR/setup.pl -d=13 -phi=3 -z=2 -g=12 -p=rmhd -eos=gamma -nf=0 -arch=gfortran10 -coord=cart

echo "grid set up complete. Now we compile"

make

echo "compilation succesful (maybe)"

mpiexec -n 4 ./bhac -i amrvac.par >output/out


echo "Compilation and run complete"

