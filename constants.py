import numpy as np
import transformations as trsfrm

# Radius of Earth in meters
Re = 6.371e6

# Magnetic moment for Earth dipole
M = 8.22e22

# Permeability of free space
u0 = 1.25663706212e-6

# Speed of light
c = 3e8

# Mass of proton in kg
M_p = 1.6726219e-27

# Mass of electron in kg
M_e = 9.10938356e-31

# Charge on electron in Coulombs
C_e = -1.60218e-19

# B scale factor for dimensionless variable
Bc = u0*M/(4*np.pi*Re**3)

# Dimensionless magnetic field at the top of the atmosphere (altitude = 100 km)
za = (Re + 1e5)/Re
Ba = np.linalg.norm(trsfrm.B(0, 0, za))