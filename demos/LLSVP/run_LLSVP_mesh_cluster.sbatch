#!/bin/bash
#SBATCH -J planet_mesh
#SBATCH -o planet_mesh_%j.txt
#SBATCH -e planet_mesh_%j.err
#SBATCH -p skx-dev
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --export=ALL
#SBATCH --time=20:00
#SBATCH -A TG-EAR170019
#SBATCH --mail-user=hrmd@mit.edu
#SBATCH --mail-type=all
./run_LLSVP_mesh.bash
