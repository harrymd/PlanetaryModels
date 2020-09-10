#!/bin/bash
cp -r /work/06414/tg857131/PlanetaryModels/ /scratch/06414/tg857131/
cd /scratch/06414/tg857131/PlanetaryModels/demos/LLSVP
mkdir -p /scratch/06414/tg857131/PlanetaryModels/output
mkdir -p /scratch/06414/tg857131/PlanetaryModels/output/LLSVP
mkdir -p /scratch/06414/tg857131/PlanetaryModels/output/unit_spheroids
module add matlab
matlab -nojvm -r "try; LLSVP_mesh_spheroidal(); catch e; fprintf(1, e.message); exit; end; exit"
