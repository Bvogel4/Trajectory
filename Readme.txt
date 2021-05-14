This program simulates a particle's motion in magnetic field, currently using a simple dipole.
2 different options for numerical integration
	Runga kutta 45, is generally faster but less accurate
	boris method, can be computed to arbitary accuracy, usually takes longer to compute. Inaccurate near low pitch angles <5deg (more than runga kutta)
Initial condition and properties are measured in a more convenient coordinate system while the integration is performed in cartesian.
An example of the motions of the plots at various pitch angles is included in the 3d plots directory. 
By default, generates trajectory at various pitch angles and plots various properties, such as error in energy, and energy in different directions along the magnetic field.
new 3d plots can be made by setting _3dplots = True when calling particle_demo. 3d animation is memory and compute intensive.
Exports a vtk for import into Paraview. 