Planetary Model Builder 
================================================================

This repository provides MATLAB scripts to build a planet model on a tetrahedal mesh,
as well as its reference gravitational field. 

How to make it work for you? 
----------------------------------------------------------------
If you only need to build a tetrahedal mesh, please compile tetgen.  
For a linux machine, just go to ./packages/tetgen1.5.0 and type make.   
If you need to compute the reference gravity, please go to packages/fmmlib3d-1.2, and compile it.  
If you use a linux machine, you may follow the README.md under ./packages

This repository uses four different packages:  
1. TetGen by Hang Si  
2. distmesh by Per-Olof Persson  
3. fmmlib3d from Leslie Greengard and Zydrunas Gimbutas  
4. a MATLAB script for vtk users by Shawn Walker

I put most of the original materials under different folders 
with some minor changes for this application.

Build your planetary models
-----------------------------------------------------------------

