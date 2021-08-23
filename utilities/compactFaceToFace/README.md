
This needs to be linked against ESI or Foundation openfoam, like:

	(source /opt/openfoam9/etc/bashrc && wmake)

Well, after changing /opt/openfoam9 according to your openfoam directory.

For installation, copy the compiled compactFaceToFace to a system PATH, 
e.g. to ~/.local/bin or ideally to /opt/openfoam9/platforms/linux64GccDPInt32Opt/bin
