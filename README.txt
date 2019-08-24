Allen Institute - Synaptic Physiology Department

This is an algorithm for determining synapse locations of neighboring human corticol neurons using 3D morphological reconstructions and
euclidean and segmental distance. This program takes in a folder of SWC files as input and returns a list of the potential synapse 
locations and a visualization of the synapses viewable in Vaa3D. Machine learning was also used to inform post-synaptic dendritic radius  
information pertinent to synapse formation.

Using the Synaptic Distance Calculator:

1) Make sure Rad_ML and Syn_Dist python files are in the same directory (as well as the Test_Files Folder)
2) Select 2 or more SWC files from the same cluster that have been (resampled at 1 pixel, have a soma estimation of 60 um, and sorted) and place them in a folder
3) run the Syn_Dist program and when prompted, select the folder with the SWC files of interest
4) after the program has run, type in "output" to receive a list of the potential synapse locations
5) run successive functions at the bottom of Syn_Dist for subsequent analysis + marker visualization 
