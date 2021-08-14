###  foamx4m
This is a minified and intermediate version of foam-extend code derived
from foam-extend-4.1, with few additions from OpenFOAM-1712 (for reading VTK files....), 
prepared by Ali Raeini for use in poreFoam codes.

* See the file Changelog.current for the list of changes.

### Licence: 
See the file COPYING for the GNU-GPLv3 copyright notice.
See also the files ExtendProjectPreamble and ListOfContributors.
The provided code is not endorsed by foam-extend developers, nor by OpenFOAM. 

### Compilling and installation
Run ./Allwmake in the top level directory after running "source ./etc/bashrc".
Once everything compiled run ./Allclean to delete intermediate files. 
An optional install script, Allinstall is also provided, but it should
be edited before running, manually setting and creating the installation paths.
The compiled libraries will have different names compared to that of OpenFOAM, 
or foam-extend, so they can be installed simultanuously, but only after all 
compilations are finished.


### notes:

For decomposepar utility to be compiled, you need to install 
libscotch-dev in your system.

Minifoam libraries do not have conflict with other versions of 
openfoam, but the executables  can be duplicate.  So if you want to 
`source $msRoot/src/script/bashrc`  while using another openfoam version, 
delete whatever you don't need from $msRoot/bin/porefoam-ext-3/ 
folder. $msRoot here stands for the upper folder where top level 
Makefile and src/ folder exist.


