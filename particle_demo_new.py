# -*- coding: utf-8 -*-
"""
Created on Sun Aug 22 17:21:12 2021

@author: blake
"""

#THis file perfroms all the transformations from the transformatons module
#computes the arrays using the Integrators module
#plots using the plots module
from datetime import datetime
startTime = datetime.now()

import numpy as np
import Transformations as trsfrm
import Integrator
import plots
import math
Re = 6.371e6
#Re = 1 # radius of earth in m
M = 8.22e22
#magnetic moment M = 8.22x10^22 \A m^2 for earth
u0 = 1.25663706212e-6

M_p = 1.6726219e-27 # kg
M_e = 9.10938356e-31
C_e = -1.60218e-19 #C


def particle_sim(L_shell = 2,
    pitchangle = 90,
    latitude = 0, #all angles in degrees for input
    longitude = 0,
    phase = 90,
    Kinetic_energy = 1e8,# in eV
    #default model = alpha particle
    mass = 4*M_p,   #in Kg
    charge = 2*-C_e,# in coulumbs
    #tp =1500#time in tc ~ 50 error gets high for rk45, boris is much more accurate
    t = .1,# time in seconds
    #1 tp corresponds to the time to go 1 radian around the gyro radius at L = 1
    accuracy = 1e2, # inverse time step in dimensionless form
    sampling = 5, # points per tc 
    method = 'rk45', # valid choices are 'boris','rk45' and 'euler'
    _2dplots = True):
    
    mname = (mass/M_p)
    qname = int((charge/C_e))
    filename = 'pitch_{}_m_{:.2}amu_q_{:}e-'.format(pitchangle,mname,qname)

    # internally all angles are in radians
    latitude = math.radians(latitude) 
    longitude = math.radians(longitude) 
    phase = math.radians(phase) 
    pitch = math.radians(pitchangle)
    #covert energy to joules
    Kinetic_energy = Kinetic_energy * 1.602176565e-19
    
    dT = 1/accuracy
    #convert initial conditons to cartesian
    x0,y0,z0, vx0,vy0,vz0 = trsfrm.ctd2car(pitch, phase, Kinetic_energy,\
                                           L_shell, latitude, longitude, mass, Re)
    S0 = np.array([x0, y0, z0, vx0, vy0, vz0])
    #convert into new dimensionless units
    #need constants first
    #re is at the top
    Bc = u0*M/(4*np.pi*Re**3)
    tc = abs((mass/(charge*Bc)))
    #in nd form can't have negative tc
    #qsign passes the sign of the charge
    #simply invert B if charge is negative.
    qsign = np.sign(charge)
    
    #convert to nd positions and velocities
    S0 = S0/Re
    S0[3:6] = S0[3:6] * tc
    #convert time into nd
    tp = t/(tc)
    #choose integrator based of input
    if method == 'rk45':
        xline,yline,zline,Vx,Vy,Vz,T = Integrator.rk45_nd(dT, tp, S0,qsign)        
    elif method =='euler':        
        xline,yline,zline,Vx,Vy,Vz,T = Integrator.euler_cromer(dT, tp, S0,qsign)
    elif method == 'boris':        
        xline,yline,zline,Vx,Vy,Vz,T = Integrator.boris(dT,sampling, S0, tp,qsign)
    else:
        print('invalid method')

    print('integrator done at time ', datetime.now() - startTime)

    if _2dplots == True:
        Bc = u0*M/(4*np.pi*Re**3)
        tc = abs(mass/(charge*Bc))
        t = T*tc
        plots.plotter(xline,yline,zline,Vx,Vy,Vz,filename, T,t,method,pitchangle,mass,Re,'s')
        print('plotter done at time ', datetime.now() - startTime)

#iterate over 2 protons, and 1 alpha particle
#electrons may need to be done alone, very memory intensive
#information on particle perdiods can be found in period_info.pdf
m = [M_p,M_p,M_e,4*M_p,4*M_p]
q = [-C_e,-C_e,C_e,-2*C_e,-2*C_e]
T = [20,20,10,30,30]
acc = [1e1,1e3,1e1,1e3,1e1]
pitch = [90,60,90,60,30]
Lshell = [2,2,8,2,2]
sample = [5,5,1,5,5]
inte = ['boris','boris','boris','boris','rk45']
#particle_sim()

#shows 5 cases showing different trajectores, particles, and integrators

for i in range(len(m)):
    particle_sim(pitchangle= pitch[i],mass = m[i], \
        charge = q[i], t = T[i],accuracy = acc[i],\
        L_shell = Lshell[i],sampling = sample[i])
