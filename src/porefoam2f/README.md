# porefoam2f module development version


 Porefoam2f code simulates incompressible two-phase flow on 3D images of porous media using Openfoam finite-volume library.



## Technical note:

### This code is experimental and needs further development:  

The interface curvature algorithm in Shams et al (2018) has lead to significant improvement in discretization of capillary forces and hence capillary pressure computation and overall accuracy of the code.  However, for complex flow domains filtering is still necessary for reducing spurious currents specially near low quality regions of unstructured meshes. Therefore, the code by Shams et al (2018) was merged with the filtering algorithm proposed in Raeini et al (2012), after major revisions.  However, when using filtering parallel to the interface, the contact line and the interface motion tends to show an artificial stick-slip behaviour. We can not use too much filtering of capillary fluxes perpendicular to the interface (which tend to alleviate the stick-slip problem), as that can lead to suppression of trapping.   In Raeini et al (2018) filtering is used as the stick-slip happens only when interfaces move, which doest not contribute as much in viscous dissipation which is what leads to flow resistance. Filtering did not severely affect the computation of capillary pressure.  If you are interested in the motion of interfaces, its is better to readjust the settings compared to default values used in the scripts here as an artificial stick-slip can affect the accuracy of interface motion.  Note stick-slip behaviour is also seen in contact line motion in nature, but for a different (physical) reason. 

If you are interested in further developing the code, please see the contact details below and get in touch.

# Instructions

Simulations can be launched using a script called AllRunImageTwoPhase. 
See [porefoam_twoPhaseFow.md](../../doc/porefoam_twoPhaseFow.md) for more details.

You can follow the series of commands run by `make test` in a copy of folder src/porefoam2f/test2f to run a simulation.
The sample simulation run with the command `make test` shows a significant level of stick-slip behaviour. 
This test case shall be updated to minimize these artefacts.

<!-- TODO: add another test case based on the star-shaped geometry used in  Raeini et al (2014). -->

<!-- TODO: porefoam2f documentation is incomplete and probably out of date. -->



### References

This code is an experimental merge of the algorithms published in he following papers:

 - M Shams, A Q Raeini, M J Blunt, B Bijeljic, “A numerical model of two-phase flow at the micro-scale using the volume-of-fluid method”, J Comp. Phys. 357:159–82 (2018) https://doi.org/10.1016/j.jcp.2017.12.027

 - A Q Raeini, M J Blunt, and B Bijeljic, “Modelling two-phase flow in porous media at the pore scale using the volume-of-fluid method”,  J Comp. Phys. 231:5653–68 (2012) https://doi.org/10.1016/j.jcp.2012.04.011


This code has been used in the following works:

 - M Shams, K Singh, B Bijeljic, M J Blunt, “Direct Numerical Simulation of Pore-Scale Trapping Events During Capillary-Dominated Two-Phase Flow in Porous Media”, Transp Porous Med (2021). https://doi.org/10.1007/s11242-021-01619-w

 - A Q Raeini, B Bijeljic, M J Blunt, “Generalized network modelling: Capillary-dominated two-phase flow”, Phys. Rev. E,  97(2):023308, (2018). https://doi.org/10.1103/physreve.97.023308
 
For other relavant publications, see:     
https://www.imperial.ac.uk/earth-science/research/research-groups/pore-scale-modelling/publications/

### Credit
  
 - Ali Q. Raeini -- developer of filtered surface force (FSF) formulation

 - Mosayeb Shams -- developer of contour-level surface force (CLSF) formulation

 - Prof. Stephane Zaleski:  supervision -- tracking interface and interfacial force computation.

 - Prof. Stephen Neethling: supervision -- sharpening of interface transition to reduce spurious currents.

 - Dr. Branko Bijeljic: PhD supervisor.

 - Prof. Martin J Blunt: PhD supervisor.


### Contacts

For more information and contacts details, see:  
https://www.imperial.ac.uk/earth-science/research/research-groups/pore-scale-modelling


### License

[GPLv3](https://www.gnu.org/licenses/gpl-3.0.txt)

<!-- TODO: add separate CREDIT.md and CHANGELOG.md -->
