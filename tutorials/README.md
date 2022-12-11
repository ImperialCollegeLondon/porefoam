# tutorials

The bash scripts [run_cfmesh.sh](run_cfmesh.sh)  or
[run_vxlmesh.sh](run_vxlmesh.sh)  runs sample simulations on a ketton
sample, which is downloaded from digital rock portal. The This is to
tests the two available mesh grneration options and the culprits of
each approach.   You need to copy them and run them from outside the
src source code folder.

[run_cfmesh.sh](run_cfmesh.sh) uses cfMesh to generate cell layers near
solid walls provided that the mesh is at least a factor of two finer
than the roughnesses of the image.   Otherwise cfMesh can crash. The
generated mesh has a high non-orthogonality but the two-phase flow
solver may be able handle that (further testing needed).


[run_vxlmesh.sh](run_vxlmesh.sh) which uses The voxelToFoam(Par)
togather with smoothMesh can in theory produce do a similar job as
snappyHexMesh (both unable to produce cell layers). Note that using
smoothMesh does not seem to improve the quality of simulations, but
this needs further testing. You can make thvoxelToFoam and
snappyHexMesh approach do not have the limitations of cfMesh but can
lead to imperfections with modelling contact line motion.

Alternatively you can edit the AllRunImageTwoPhase script invoked by
run_vxlmesh.sh and uncomment the section which uses snappyHexMesh. This
requires installation of a full openfoam version and of course further
testing.  Another area which needs further testing is the smoothMesh utility
whose mesh smoothing parameters are not roboust and need adjustement in
a case by case basis.
