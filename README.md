# Overview

This program simulates a particle's non-relativistic motion in dipole magnetic field.

Three different options are available for numerical integration: [Euler-Cromer](https://aapt.scitation.org/doi/10.1119/1.12478), [Boris](https://books.google.com/books?id=S2lqgDTm6a4C&q=Borris#v=onepage&q=Boris&f=false), and [Runge-Kutta-Fehlberg](https://ntrs.nasa.gov/api/citations/19700031412/downloads/19700031412.pdf).

Euler-Cromer integration is faster, but less accurate than both Boris and RKF.

`rk45` takes longer for similar accuracy, except when adaptive time step is especially useful such as at high latitudes.

The Boris method can be computed to an arbitrary time step, usually shorter for similar accuracy. Inaccurate near small pitch angles (< 5°) (more than Runga Kutta) due to no adaptive steps. 

For equatorial-bound particles, the time step needed for the Boris method is quite low at `0.01` in nd time steps otherwise lower timesteps are generally needed. The Boris integration option in this code also has an option to only simulate one drift period in order to shorten computation time. 

The Boris method is tailored specifically to Lorentz force motion, as it conserves system's energy better. Boris method in general can handle electric fields, but this has not been implemented.

Initial condition and properties are measured in a more convenient centered dipole coordinate system while the integration is performed in non-dimensional cartesian coordinates.

# Use
```
git clone https://github.com/rweigel/Trajectory
cd Trajectory
python3 command
```

To to modify run, edit top of `ParticleDemo.py`.



# Demo Example
Function demo generates various verification plots of to verify accuracy and compare algorithims. This generates various plots with different parameters for protons and elections, specifically L-shell, and kinetic energy. The results of these simulations can be found in [this folder](plots) 
The function trajectory_generator generates ~800 various trajectories for import into paraview, since one the fly generation is too computationally expensive

With _2dplots = True, generates 3 plots: error vs time, cartesian coordinates vs time, and x vs y
more plots can be enabled by uncommenting them in plots&#46;py. For checking accuracy or long simualtions, the energy.svg are crucial to checking to make sure the simulation maintains the desired level of accuracy.

Exports particle path trajectories as .VTKs for import into Paraview for viewing.
# Comparison of algorithims
![Boris Method](plots/qm_95,788,534.15706204_Ke_100000.0MeV_pitch_90d_L_3Re_boris/energy.svg)
![Runge Kutta 4th order](plots/qm_95,788,534.15706204_Ke_100000.0MeV_pitch_90d_L_3Re_rk45/energy.svg)


## Appendix

A sample plot of a proton Gyrofrequency, and the corresponding verification plot
![Proton Gyro](plots/qm_95,788,534.15706204_Ke_100000.0MeV_pitch_90d_L_1.64Re_boris/x.svg)
![Gyrofrequecny and gyro period of electrons and protons](Period_info/Period_info-1.png)

