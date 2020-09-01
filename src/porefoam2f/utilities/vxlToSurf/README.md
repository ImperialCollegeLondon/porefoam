### vxlToSurf  - convert voxelImage (.mhd, .tif, .am files) to a surface mesh (.vtk/.obj file)


This folder contains parts of a library I am developing for conversion of 3D images into surface meshes representing their boundary and their processing and analysis.  

It is included in porefoam code  as it is a pre-requisite for mesh generation using snappyHexMesh and cfMeshs, invoced through the AllRunImageTwoPhase* scripts.

The code depends on libvoxel but not on other openfoam codes.
