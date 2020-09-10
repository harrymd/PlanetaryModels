# harrymd/PlanetaryModels

This repository is a fork of [Jia Shi's *PlanetaryModels*](https://github.com/js1019/PlanetaryModels) code. It contains scripts for building Earth models containing [large low-shear-velocity provinces](https://en.wikipedia.org/wiki/Large_low-shear-velocity_provinces) (LLSVPs).

If you wish to use the examples here, you will have to substitute your own file paths (for example, `/your/work/directory` instead of `/work/06414/tg857131`).

## Building PREM and LLSVP models

Move to the script directory:

```bash
cd /work/06414/tg857131/PlanetaryModels/demos/LLSVP/
```

The model is controlled by two parameters as described by `input_format.txt`:

```
tet_max_vol (km3)   Maximum volume of tetrahedron allowed.
pOrder              Order of finite elements.
```

These parameters are set in `input.txt`, for example:

```
1.0E9
1
```

Then simply

```bash
./run_LLSVP_mesh.bash
```

Or, if it is a large model, use the cluster:

```bash
sbatch run_LLSVP_mesh_cluster.sbatch
```

The input and scripts are copied to `$SCRATCH` before execution. Note that any output files already in `$SCRATCH` with the same input parameters will be overwritten. To view the output files:

 ```bash
 cd /scratch/06414/tg857131/PlanetaryModels/output
 ```
 
The new model files will be found in `LLSVP/` in directories such as `llsvp_2039.6_2.00_1`. New unit sphere files such as `unit_sphere_1.100000E+00.mat` will be found in `unit_spheres/`. To send the new output files to the portable hard drive, run:
 
```bash
rsync -rvh tg857131@stampede2.tacc.utexas.edu:/scratch/06414/tg857131/PlanetaryModels/output/ /Volumes/stoneley5TB/all/PlanetaryModels/v2
```

and to send them to the local machine, use:

```bash
rsync -rvh tg857131@stampede2.tacc.utexas.edu:/scratch/06414/tg857131/PlanetaryModels/output/ /Users/hrmd_work/Documents/research/stoneley/output/PlanetaryModels
```

You can then check how they look by opening the `.vtk` files in *Paraview*.

## Possible errors
 
You may encounter
 
```
./tetgen: cannot execute binary file
```
 
which can be fixed by recompiling *Tetgen* (simply `make clean` and `make`).