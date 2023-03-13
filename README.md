# Diffraction Post-processing Tools for Crystal Plasticity Simulations Conducted Using Matflow
This repository contains simulation post-processing python functions to calculate and plot the properties of a crystal plasticity simulation, analagous to the quantities that can be obtained from XRD analysis of the microstructure.
We aim to create a set of tools to make comparison of experimental XRD results and simulation results easily comparable for model validation and exploration of deformation mechanisms.
Using a method to mask voxels of the simulation which satisfy Bragg's condition for a given axis direction (X,Y,Z), stress, strain and slip activity results can be filtered out for specified crystallographic planes of the simulated phase.

See [here](https://lightform-group.github.io/wiki/software_and_simulation/matflow-post-processing) for more guidance on matflow results and methods to explore/analyse simulation results using matflow's python tools.

Update 2023-03-07:
The notebook has been discussed in the monthly Lightform modelling meeting and has been made publicly available to all to review/add to. If you have any ideas for additional post-processing simulation results you'd like to include please submit an issue or pull request.
