# Overview

This program simulates a charged particle's relativistic motion in dipole magnetic field.

Three different options are available for numerical integration: [Boris](https://books.google.com/books?id=S2lqgDTm6a4C&q=Borris#v=onepage&q=Boris&f=false), [Runge-Kutta-Fehlberg](https://ntrs.nasa.gov/api/citations/19700031412/downloads/19700031412.pdf), and [Euler-Cromer](https://aapt.scitation.org/doi/10.1119/1.12478)


`rk45` ues the Python [SciPy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.RK45.html) implentation.

The Boris method can be computed to an arbitrary time step and is more compuationally efficient than `rk45`. The Boris method maintains accuracy up to steep pitch angles, but does not use a true adaptice timestep as the error is predictable with changes in latitude.

For equatorial-bound particles, the time step needed for the Boris method is quite forgiving around `0.01` in non-dimeionsal time steps otherwise lower timesteps are generally needed. The Boris integration option in this code also has an feature to only simulate one drift period in order to shorten computation time. 

The Boris method is tailored specifically to Lorentz force motion, as it conserves system's energy better. Boris method in general can handle electric fields, but this has not been implemented.

Initial conditions are specified in radius, latitude, longitude, pitch angle, kinetic energy, and phase angle.

# Use

```
git clone https://github.com/Bvogel4/Trajectory
cd Trajectory
python3 
```

# Main Functions

There are three main functions that can be used:

1. `particle_sim()`, which generates the arrays containing the positions, velocites and times of a particle.
2. `save()`, which writes exports particle path trajectories as .TXT and .VTKs, the latter for import into Paraview for viewing.
3. `plot()`, which generates a list of plots including: cartesion postions and velocites, absolute error in energy, and an xy top-down view.

# Demos and Examples

1. `demo_method()` generates various verification plots to verify accuracy and compare algorithims. This generates various plots with different parameters for pitch angle. (the primary driver for differences in error).
2. `demo_bounce()` computes a full bounce for an electron and a proton in an near-equatorial orbit. 
3. `demo_drift()` computes a full drift period for an electron and a proton in an equatorial orbit. 
4. `trajectory_generator()` generates ~800 various trajectories for import into paraview, since on the fly generation is too computationally expensive. 

The results of these demos will be be found in [this folder](output).

# Comparison of algorithims

To compare methods a metric that measure the error in energy per unit time was used, called compute efficeincy. Comparing this value for differnt methods and pitch agnles serces to quantify the quality of the algorithims. (lower value is better)

| Pitch angle | Boris .01 | Boris .0001 | rk45   |
|-------------|-----------|-------------|--------|
| 90          | 6.3e-14   | 3.4e-15     | 3.5e-9 |
| 10          | 2.1       | 2.9e-9      | 4.3e-9 |

As can be seen from the table, the Boris method has significantly higher accuracy per unit computation time than `rk45` at equatorial orbits (pitch angle of 90 degrees), and still remains competitive at low pitch angles. 
