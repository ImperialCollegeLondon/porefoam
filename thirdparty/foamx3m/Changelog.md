
************************************************************************************
This file contains brief description of the changes in current version of the code
It can also contain the changes from previous releases. The general format should be:


### The previous changes should be marked with "Release 4.0 2018/03/2018", for instance.

## The changes should be in one of the two section tags: Bugfix and Feature

* Branchname
  One or two lines describing the bugfix or feature
  Author: Contributor name; Merge: Code maintainer name (if different from Author)

***********************************************************************************


## TODO
*  Problem to resolve: all the /*- ... -*/ comments in the middle of source 
   codes are deleted, find a way to put them back. 
*  the utilities provided don't work with poerFoam

## Bugfix
* Rolling back from foam-extend-4.0 to 3.1: 
  The finite-volume library and most of parallelization routines are replaced with 
  those from foam-extend-3.1 Branch.  This was done to avoid a crash in openfoam related 
  to a CFD-Online froum post title "OpenFOAM-2.1 instability problems?" on May 25, 2012
  Author: Ali Raeini

## Feature
* Compatibility with OpenFOAM
  added ref() functions to tmp and GeometricField classes.
  argList class changed to allow some of the openfoam application-specific 
  options and arguments
  Author: Ali Raeini

* Minified foam-extend: 
  Removed third party, most of the applications/utilities, all the solvers.
  The results is a minimal foam-extend library usefull for those who need
  the standard parts of the openfoam library. The shell commands run to help 
  the merge are recorded in the file foam_clean_commands.
  Author: Ali Raeini



### Release 4.0 2018/03/2018

## Feature

## Bugfix
* updatePreconfiguredBoundaryTuts
Author: Hrvoje Jasak; Merge: Dominik Christ
Updated tutorials which use mixing plane and had broken preconfigured boundary
files

* ReleaseNotesInstallationVersionNumber
Author: Hrvoje Jasak; Merge: Dominik Christ
Bumped version number in ReleaseNotes and updated links to SourceForce pages

* Intel14Port
Author: Hrvoje Jasak; Merge: Dominik Christ
Updated wmake/rules for compilation with ICC 14

* processorPointPatch
Author: Zeljko Tukovic; Merge: Dominik Christ
Correct interpolation of point data on processor boundaries after parallel
topological changes


