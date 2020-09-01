### vxlToSurf  - convert voxelImage (.mhd, .tif, .am files) to a surface mesh (.vtk/.obj file)


This folder contains parts of a surface library that I am developing for conversion of 3D images into surface meshes representing their boundary and their processing and analysis.  

It is included in porefoam code  as it is a pre-requisite for mesh generation using snappyHexMesh and cfMeshs, invoked through the AllRunImageTwoPhase* scripts.  


All the codes here compile into a single executable with the same name as the directory name.  The codes depend on libvoxel but not on other openfoam codes.
They replace a previous code named voxelToSurface which was written using openfoam libraries and hence too difficult to integrate with my other C++ codes.
