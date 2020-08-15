# porefoam2f module development version


 Porefoam2f code simulates incompressible two-phase flow on 3D images of porous media using Openfoam finite-volume library.



## Technical note:

### This code is experimental and needs further development:  

The interface curvature algorithm in Shams et al (2018) has lead to significant improvement in discretization of capillary forces and hence capillary pressure computation and overall accuracy of the code.  However, for complex flow domains I realized filtering is still necessary, that is why I merged the code by Shams et al (2018) with the filtering algorithm proposed in Raeini et al (2012), well after major revisions.  The Main challenge with using filtering is that when using filtering parallel to the interface, the contact line and the interface, tend to have a artificial high stick-slip behaviour. We can not use too much filtering of capillary fluxes perpendicular to the interface (which tend to alleviate the stick-slip problem) neither, as that can lead to suppression of trapping.   In Raeini et al (2018) I use filtering as the stick-slip happens only when interfaces move, which doest not contribute as much in viscous dissipation which is what leads to flow resistance. Filtering did not severely affect the computation of capillary pressure.  If you are interested n the motion of interfaces, you probably better to use different settings than recommended here as, obviously, an artificial stick-slip makes it difficult to predict interface motion accurately.  Note stick-slip behaviour is also seen in contact line motion in nature, but for a different (physical) reason. 

We aim to improve boundary conditions and interface tracking algorithm, which has been put on hold due to lack of free time / dedicated funding.  See below for contact details if you want to help with your time or funding.

# Instructions

You can follow the series of commands run by `make tst` in a copy of folder src/porefoam2f/test2f to run a simulation.
The sample simulation run with the command `make tst` shows a significant level of stick-slip behaviour. In future, I might update the filtration parameters etc to produce a better test case.  

TODO: add another test case based on the star-shaped geometry used in  Raeini et al (2014). 

TODO: porefoam2f documentation is incomplete and probably out of date.



### References

This code is an experimental merge of the algorithms published in he following papers:

 - M Shams, A Q Raeini, M J Blunt, B Bijeljic, “A numerical model of two-phase flow at the micro-scale using the volume-of-fluid method", J Comp. Phys. 357:159–82 (2018)

 - A Q Raeini, M J Blunt, and B Bijeljic, “Modelling two-phase flow in porous media at the pore scale using the volume-of-fluid method”,  J Comp. Phys. 231:5653–68 (2012)


The code will not produce identical results to neither of the above papers though.  The results in the following paper are produced using this code, version from 2018, which has not been changed much since then until the first Github release of this code (2020):

 - A Q Raeini, B Bijeljic, M J Blunt, “Generalized network modelling: Capillary-dominated two-phase flow”, Phys. Rev. E,  97(2):023308, (2018)
 
For other relavant publications, see:     
https://www.imperial.ac.uk/earth-science/research/research-groups/perm/research/pore-scale-modelling/publications/

### Contacts

For any queries, contact:     
 - Ali Q. Raeini, email: a.q.raeini@imperial.ac.uk
 - Mosayeb Shams, email: m.shams14@imperial.ac.uk

For more information and contacts details, see:  
http://www.imperial.ac.uk/earth-science/research/research-groups/perm/research/pore-scale-modelling  

### License

[GPLv3](https://www.gnu.org/licenses/gpl-3.0.txt)
