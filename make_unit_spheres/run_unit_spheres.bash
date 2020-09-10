#!/bin/bash

#SBATCH -J unit_sphere
#SBATCH -o unit_sphere_%j.txt
#SBATCH -e unit_sphere_%j.err
#SBATCH --time=2:00:00
#SBATCH -A TG-EAR170019
#SBATCH --mail-user=hrmd@mit.edu

matlab -nojvm -r "tet_max_vol=8.0E10; radius=6371.0; out_dir = './output';  distmesh_path = '../distmesh'; try; run_unit_sphere(tet_max_vol, out_dir, distmesh_path, radius); catch; exit; end; exit"