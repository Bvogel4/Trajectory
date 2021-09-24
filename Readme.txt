This program simulates a particle's non-relativistic motion in earth's magnetic field, modeled as a simple dipole. Can easily be modified for any dipole. 
3 different options for numerical integration:
Euler Cromer is faster, but less accurate than both Boris and RK45.

Runga Kutta 4th order (rk45) takes longer for similar accuracy, except when adaptive time step is especially useful such as at high latitudes.

Boris method can be computed to an arbitrary time step, usually shorter for similar accuracy. Inaccurate near low pitch angles <5deg (more than Runga Kutta) due to no adaptive steps.
For equatorial bound particles in Boris method, the time step need is quite low at 1e-1 in nd time steps otherwise lower timesteps are generally needed. The Boris method also has a feature to only simulate 1 drift period to shorten computation time. 
Boris method is tailored specifically to Lorentz force motion, as it conserves system's energy better. Boris method in general can handle E fields, but this functionality is currently omitted.

Initial condition and properties are measured in a more convenient centered dipole coordinate system while the integration is performed in non-dimensional cartesian coordinates.
#An example of the motions of the plots at various pitch angles is included in the 3d plots directory. 
By default, generates trajectory of electrons, proton, and alpha particles at a few varying pitch angles. Note that electrons are particularly compute intensive due to the high charge to mass ratio
with _2dplots = True, generates 3 plots: error vs time, cartesian coordinates vs time, and x vs y
more plots can be enabled by uncommenting them in plots.py

Exports particle path trajectories as .VTKs for import into Paraview for viewing.
