# -*- coding: utf-8 -*-
"""
Created on Fri Feb 25 12:49:01 2022

@author: blake
"""
import numpy as np
import transformations as trsfrm

# constants
Re = 6.371e6
# Re = 1 # radius of earth in m
M = 8.22e22
# magnetic moment M = 8.22x10^22 \A m^2 for earth
u0 = 1.25663706212e-6
c = 3e8
# critical B value
Bc = u0*M/(4*np.pi*Re**3)

M_p = 1.6726219e-27  # kg
M_e = 9.10938356e-31
C_e = -1.60218e-19  # C
# magnetic field at the top of the atmopshere in nd terms
za = (Re + 1e5)/Re
Ba = np.linalg.norm(trsfrm.B(0, 0, za))