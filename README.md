

 ----------------------------------------------------------------    
 
## This repository is released in the hope that it will be useful for experts in direct two-phase flow simulation and is not intended for everyday use.  See  [src/porefoam2f](src/porefoam2f) for more specific details.

See also README files for other modules  which are located in their own directories. 

 ----------------------------------------------------------------

## General notes

### Compiling

To compile, open a terminal in the upper most directory and run:    

 `make -j`

To test the compilation, run:    

 `make test`

Once everything compiled successfully, to clean the temporary files, type:

 `make clean`

The above command can be run inside most of the subfolders, wherever a 
makefile or Makefile is present.  The libraries, those with a `makefile`,
should be compiled before the apps that contain `Makefile`s.

Compilation requires gnu (Linux) make, cmake, a c++ compiler with -std=c++11
support and an MPI. The compilation is tested using g++ (version 5+) (default)
and using intel-2018 compilers.


### Tests and demos
To test the codes type:

 `make test`

This should copy a series of input files/scripts in a `test` folder and 
run a series of relatively quick test cases (see README.md files in 
subdirectories).  

### Coding conventions

GNU Makefile scripts are used primarily for code compilation and 
running quick tests by code developers.

Automatic tests are written using input files for C++ codes, new C++ 
executables testing internals of the codes and Python scripts. All 
these can be run using `make test` command, which uses 
script/testApp bash script.

All scripts, either for testing or production, which need mathematical 
calculations or plotting and are not performance critical are developed 
using Python. We use Python 3, shich should be available as python3 command.
In Ubuntu (18.04) this can be installed by typing in a terminal:    
 `sudo apt install python3`

All computations which typically take more than few minutes are 
developed using C++.  In some variants of the code, where interaction 
with other languages (R/Shell/Python) are needed, the C++ code shall 
make use of SiR class to achieve this, where possible, using the `_sh_` 
command to launch other scripts/executables.


Bash/Shell scripts are used to launch run Makefile and Python scripts, 
either for testing or in production to simplify the run of openfoam 
solvers.  



### Contact and References ###

For contacts and references please see:    
http://www.imperial.ac.uk/earth-science/research/research-groups/perm/research/pore-scale-modelling    
or contact:    
 Ali Q. Raeini, email: a.q.raeini@imperial.ac.uk     
 Mosayeb Shams, email: m.shams14@imperial.ac.uk     

More details are given in the src/doc directory.

