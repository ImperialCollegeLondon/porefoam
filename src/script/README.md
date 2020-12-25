# Make/installation scripts

This module contains bash and makefile script for helping with Cmake, 
GNU- make and Bash script required to compile and run different 
**lib**rary and **app**lications developed by Ali Q. Raeini (2020).  

The `initbash` is an independent bash script containing utility macros, 
which together with `bashrc`  installation script, is required for 
running workflows developed as bash scripts.

`msrc.py` is needed for running workflows written in python3 scripts.

---------

## Aims 

These scripts has been released as a separate module, to help 
the code developers with script re-use, separation of content from 
logic, and simplification of workflow for code compilation and release. 
Hopefully it can help others as well.

---------

##  Usage and caution

The script here use recursive make by running the AllMake and AllClean 
scripts. The script change, add and delete and delete files on your 
computer: the root directory where files are changed, generated or 
deleted is called msRoot which, by default, points to two directories 
upper to the location of these script themselves.  If you ever think of 
using these scripts for building your applications, make sure these are 
wrapped inside two (sub-)sub-folders, dedicated for code development, 
with the important (source) files regularly backed up.  Here is how 
the directory structure should look like:


- `apps/ -------------------- msRoot directory`

    - `src/ ----------------- -- source codes`
        * `script/ ---------- -- ** this module, build & install`
        * `include/ --------- -- ** C++ utility codes`
        * `bench/ ----------- -- ** test/example data`
        * `...`
        * `...`
        
    - `thirdparty/ ---------- -- others' source codes`
        * `foamx3m ---------- -- ** a minified openfoam `
        * `svplot ---------- -- ** a modified former svg_plot`
        * `zlib`
        * `libtiff`
        * `...`

    - `bench/ --------------- -- large input files, too large for src/`
    - `test/ ---------------- -- test folder (auto copied from bench/)`
    - `bin/ ----------------- -- executables folder (auto crea/deleted)`
    - `lib/ ----------------- -- library files (auto crea/deleted)`
    - `share/ --------------- -- configuration files (auto crea/deleted)`
    - `build/ --------------- -- temporary files (auto crea/deleted)`
    - `...`


As shown above, this `script` folder is typically located in `src` 
(which contains regularly changed source codes as subdirectories). The 
less frequently changed source code are placed a directory called 
`thirdparty`.  A bench folder is also sometimes included holding 
temporary/client data, as shown above.  All `.git` directories are kept 
outside the msRoot directory and are used to test and release varius 
modules of the code, which end up in the msRoot directory in the 
modules published using git.  The `src` folder is regularly backed up,
it is recommended to make regular backups of `thirdparty` and 
`bench` folders as well.


---------

## TODOs

Convert all `All*` bash script into makefile script (?).

Make use of cmake more, for cross-platformt portability. Explore possibly 
of having two make versions: gnu makefile script and a Cmakelists.txt.

---------
