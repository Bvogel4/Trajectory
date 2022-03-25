# Overview

This program simulates a charged particle's relativistic motion in dipole magnetic field.

Three different options are available for numerical integration: [Boris](https://books.google.com/books?id=S2lqgDTm6a4C&q=Borris#v=onepage&q=Boris&f=false), [Runge-Kutta-Fehlberg](https://ntrs.nasa.gov/api/citations/19700031412/downloads/19700031412.pdf), and [Euler-Cromer](https://aapt.scitation.org/doi/10.1119/1.12478)

`rk45` ues the Python [SciPy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.RK45.html) implentation.

The Boris method is tailored specifically to Lorentz force motion, as it conserves system's energy better. THe Boris method in general can handle electric fields, but this has not been implemented.

The Boris Method implementation is generally more effiecient than 'rk45', and uses a default, or user set non-dimeionsal (ND) time-step. The implemenation of the boris method does not feature a true adaptive time step, but instead a technique that decreases the ND time-step as the lattitude approaches 90 or -90. This is the main factor that contributest to different error sizes.

For equatorial-bound particles, the time step needed for the Boris is around `0.01` in ND time steps otherwise lower timesteps are generally needed. The Boris implementation is designed to only simulate one drift period in order to shorten computation time. 

Initial conditions are specified in radius, latitude, longitude, pitch angle, kinetic energy, and phase angle.

# Use

```
git clone https://github.com/Bvogel4/Trajectory
cd Trajectory
python3 demo_method.py 
python3 demo_bounce.py 
python3 demo_drift.py 
python3 trajectory_generator.py 
```

# Main Functions

There are three main functions that can be used:

1. `particle_sim()`, which generates the arrays containing the positions [RE], velocites[RE/s] and times[s] of a particle.
2. `save()`, which writes exports particle path trajectories as .TXT and .VTKs, the latter for import into Paraview for viewing.
3. `plot()`, which generates a list of plots including: cartesian postions and velocites, absolute error in energy, and an xy top-down view.

# Demos and Examples

1. `demo_method()` generates various verification plots to verify accuracy and compare algorithims. This generates various plots with different parameters for 2 pitch angles. (the primary driver for differences in error).
2. `demo_bounce()` computes a full bounce for an electron and a proton in an near-equatorial orbit. 
3. `demo_drift()` computes a full drift period for an electron and a proton in an equatorial orbit. 
4. `trajectory_generator()` generates ~800 various trajectories for import into paraview, since on the fly generation is too computationally expensive. 

The results of these demos will be be found in [this folder](output).

# Comparison of algorithims and Verification/Validation

To compare methods a metric that measure the error in energy per unit time was used, called compute efficeincy. Comparing this value for differnt methods and pitch agnles serces to quantify the quality of the algorithims. (lower value is better)

See also [demo_method.py](https://github.com/Bvogel4/Trajectory/blob/main/demo_method.py).

| Pitch angle | Boris .01 | Boris .0001 | rk45   |
|-------------|-----------|-------------|--------|
| 90          | 6.3e-14   | 3.4e-15     | 3.5e-9 |
| 10          | 2.1       | 2.9e-9      | 4.3e-9 |

As can be seen from the table, the Boris method has significantly higher accuracy per unit computation time than `rk45` at equatorial orbits (pitch angle of 90 degrees), and still remains competitive at low pitch angles. 
