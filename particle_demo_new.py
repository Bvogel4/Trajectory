# -*- coding: utf-8 -*-
"""
Created on Sun Aug 22 17:21:12 2021

@author: blake
"""
from datetime import datetime
startTime = datetime.now()
import numpy as np
import Transformations as trsfrm
import Integrator as Integrator
import plots
import math
Re = 6.37e6
#Re = 1 # radius of earth in m
M = 8e22 #magnetic moment M = 8x10^22 \A m^2 for earth
u0 = 1.25663706212e-6
#for j in [90,60,45,30,0]:
    # i could add multithreading to this step for multiple generations
for pitchangle in [90]: # generates plots with different intial pitch angles
    
    L_shell = 1
    latitude = 0 #all angles in degrees
    longitude = 0
    phase = 90
    Kinetic_energy = 1e6 # in eV/kg
    #model of an alpha particle
    mass = 6.68e-27  #in Kg
    charge = 3.2e-19  # in coulumbs
    tp = 1#time in tc ~ 50 error gets high for rk45, boris is much more accurate
    #1 tp corresponds to the time to go 1 radian around the gyro radius at L = 1
    accuracy = 1e3 #  just inverse time step
    sampling = 5 # points per tc 
    method = 'boris'
    _2dplots = True
    
    filename = 'pitch{0}'.format(pitchangle)
    
    
    
    # internally all angles should be in radians
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
    #re is at he top
    Bc = u0*M/(4*np.pi*Re**3)
    tc = mass/(charge*Bc)
    #convert to nd positions and velocities
    S0 = S0/Re
    S0[3:6] = S0[3:6] * tc
    
    #non dimenionlize r. t is already nd
    #S0 = np.array([1,0.,0.,0.,.01,.01]) # temp override of initial conditions
    
    
    
    if method == 'rk45':
        
        #tc = mass/(charge*Bc)
        
        #print(Bc,dT)
        #nd initial conditions
        #S0[0:3] = S0[0:3]/Re
        
        #print(S0)
        xline,yline,zline,Vx,Vy,Vz,T = Integrator.rk45_nd(dT, tp, S0)
        #convert back into dimensional coords?
        
    elif method =='euler':
        
        xline,yline,zline,Vx,Vy,Vz,T = Integrator.euler_cromer(dT, tp, S0)
        #print(T)
        #print(xline,yline)
    elif method == 'boris':
        
        xline,yline,zline,Vx,Vy,Vz,T = Integrator.boris(dT,sampling, S0, tp)
    else:
        print('invalid method')

    print('integrator done at time ', datetime.now() - startTime)

    if _2dplots == True:
        Bc = u0*M/(4*np.pi*Re**3)
        tc = mass/(charge*Bc)
        t = T#*tc
        plots.plotter(xline,yline,zline,Vx,Vy,Vz,filename, T,method,pitchangle,mass,Re)
        print('plotter done at time ', datetime.now() - startTime)
        
        
            #should do conersions after function? no I should do it before and after calling 

    #particle_demo(pitch =j,_3dplots=False,fft = False, method = 'nd',_2dplots=(True))