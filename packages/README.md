External packages
------------------------------------------------------------------
Here are four packages that are used for model building:  
1. distmesh by Per-Olof Persson.   
2. fmmlib3d from Leslie Greengard and Zydrunas Gimbutas.    
3. A Matlab script for vtk users by Shawn Walker.  
4. TetGen by Hang Si.

### 1. distmesh

Usually no compilation required. If you receive an error message about invalid MEX files, it may help to re-compile the `.cpp` functions. For example, to recompile `trisurfupd.cpp`, use the following commands in the Matlab shell:

```matlab
mex -setup C++
mex trisurfupd.cpp
```

### 2. fmmlib3d-1.2

Fast Multipole Method. You only need to compile this if you want to compute the reference gravity.

#### 2.1 Linux

For a Linux machine with Matlab installed, please do 

```bash
make lib
make mwrap
make mex-matlab 
```

check `make test-mex-matlab` and make sure that it works.
You can then call it in the Matlab scripts.

#### 2.1 Mac

Same as above.

If you get an error 'Actual argument contains too few elements for dummy argument', use the the Fortran flag `-std=legacy`. E.g.

```bash
make FFLAGS='-std=legacy' lib
```

To compile the matlab components (`mwrap`, `test-mex-matlab`) you need to be able to call Matlab from the command line:

```bash
export PATH=$PATH:/Applications/MATLAB_R2018b.app/bin/
```

You cannot `make all` unless you have Octave installed.

If Matlab has a problem finding a `*.dylib file`, but you have the file somewhere on your system, a quick solution is to place a symbolic link to the file in the place where Matlab is looking, e.g.

```bash
ln -s /opt/local/lib/libgcc/libgfortran.3.dylib /usr/local/lib/libgfortran.3.dylib
```

#### 2.2 Other platforms

Not tested on other platforms. You may have to modify the makefiles.  

### 3. VTK scripts

No compilation required.

### 4. TetGen

```bash
cd packages/tetgen1.5.0 
make
```
