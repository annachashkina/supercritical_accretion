# supercritical_accretion
The code calculates a supercritical accretion disc structure around a magnetized neutron star. 
Basic physics is described in two papers, http://adsabs.harvard.edu/abs/2017MNRAS.470.2799C and Chashkina et al., in preparation. 
Most of the physics is given in {\tt physics.py} and {\tt dirin.py}, in the latter file, solution is found for the eigenvalue problem. 
{\tt meshes.py} calculates the grids for different parameter ranges
{\tt plots.py } performs reprocessing and visualization.
{\tt readkit.py } reads the data written in datafiles.
{\tt ssdisk.py } contains Shakura-Sunyaev's approximate analytic solutions.
