import numpy as np
import transformations as trsfrm

# Radius of Earth in meters
Re = 6.371e6

# Magnetic moment for Earth dipole
M = 8.22e22

# Permeability of free space
u0 = 1.25663706212e-6

# Speed of light in m/s
c = 299792458

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



#default parameters
parameters = {'L_shell':2,'pitch_angle' : 90,'mass' : M_p, 
              'charge' : -C_e,'Kinetic_energy':1e8,'time':1,
              'species': 'proton','latitude':0,'longitude':0,'phase':0,
              'method' :'boris', 'accuracy' : 1e3,'sampling':36,
              'loss_cone':True,'show_timing':False}

'''
L_shell=2,
pitch_angle=60,
mass=constants.M_p,     # in kg
charge=-constants.C_e,  # in Coulumbs
latitude=0,     # all angles in degrees
longitude=0,
phase=0,
Kinetic_energy=1e8, # in eV. Defaults to high energy for 
                    # short drift period

t=1e1,          # time in seconds

method='boris', # valid choices are 'boris','rk45'
                # and 'euler_cromer'

accuracy=1e3,   # inverse time step in dimensionless form
sampling=36,    # points per gyro
                # note accuracy*sampling cannot be greater than 1

losscone=True,  # True to ditch atmoshperic particles
                # False to keep them, like for demo()
show_timing=False
'''