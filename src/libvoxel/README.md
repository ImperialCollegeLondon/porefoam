##  libvoxel

This library serves as a header-only 3D image manipulation and I/O for my other 
software which work with X-ray computer tomography data.

In addition, voxelImageProcess and voxelToFoam(Par) applications are included here for convinience; the libvoxel codes (all starting with voxelImage) are independent of these two packages.


The library can read raw data in ascii (.dat) or binary (.raw) formats, in Avizo (.am) formats (only uncompressed and ByteRLE encoded data are supported).  It can also read raw.gz and .tif image formats provided that the  [libz] and [libtiff] libraries are available.

### Usage

This library is used to read 3D image files from other codes, however the standalone app `voxelImageProcess`, which solely acts as an interface to libvoxel, can be used to print help messages about the keywords supported by libvoxel:

   `voxelImageProcess -h`


### Compile instructions

#### pre-requisites:

To install necassary libraries in Ubuntu (18.08 etc.), run:    
	`sudo apt install libjpeg-dev  #used by libtiff`    
	`sudo apt install liblzma-dev  #used by libtiff`    

### Download, build and install
This library is not to be compiled (after all its a header-only template library), and a component of my other `apps`: for instance [pnextract]/[pnflow] and [porefoam].  

Nevertheless, if you want to generate the voxelImageProcess and voxelToFoam(Par) executables, run the following commands in order (tested in Linux Ubuntu 18.04)

```shell

    git clone https://github.com/aliraeini/script apps

    cd apps

    git clone https://github.com/aliraeini/libvoxel src/libvoxel

    git clone https://github.com/aliraeini/libtiff thirdparty/libtiff

    git clone https://github.com/aliraeini/zlib thirdparty/zlib

    make -f src/script/Makefile.msRoot
```

This library is not actively developed by itself, it is however kept syncronised with a not released development branch.  We will be happy to merge and keep record of any pull-requests you send us though.

###  Licence

The code and executables are provided as is, without any kind of warranty;
use at your own risk.

For further information contact me by email:   a.q.raeini@imperial.ac.uk


### References
See the [Publications on our website], also [Images on our website].

[Publications on our website]: https://www.imperial.ac.uk/earth-science/research/research-groups/pore-scale-modelling/publications/
[Images on our website]: https://www.imperial.ac.uk/earth-science/research/research-groups/pore-scale-modelling/micro-ct-images-and-networks/
[Imperial College - pore-scale consortium]: https://www.imperial.ac.uk/earth-science/research/research-groups/pore-scale-modelling
[libtiff]: https://gitlab.com/libtiff/libtiff
[porefoam]: https://github.com/aliraeini/porefoam
[pnextract]: https://github.com/aliraeini/pnextract
[pnflow]: https://github.com/aliraeini/pnflow
[libz]: https://github.com/madler/zlib
