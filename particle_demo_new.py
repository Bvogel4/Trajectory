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
for pitchangle in [45]: # generates plots with different intial pitch angles
    filename = 'pitch{0}'.format(pitchangle)
    L_shell = 4
    latitude = 0 #all angles in degrees
    longitude = 0
    phase = 90
    Kinetic_energy = 3e10 # in eV/kg
    mass = 1e-16  #in Kg
    charge = 5e-18  # in coulumbs
    tp = 2000000#time in tc ~ 50 error gets high for rk45, boris is much more accurate
    #1 tp corresponds to the time to go 1 radian around the gyro at L = 1
    accuracy = 1e2
    sampling = 1 # points per tc 6tc ~ around 
    method = 'boris'
    _2dplots = True
    # internally all angles should be in radians
    latitude = math.radians(latitude) 
    longitude = math.radians(longitude) 
    phase = math.radians(phase) 
    pitch = math.radians(pitchangle)
    
    Kinetic_energy = Kinetic_energy * 1.602176565e-19
    
    #convert initial conditons to cartesian
    
    #is ctd2car working right? I don't think it is
    #I think it's mostly working, I'll cross check with kyle to make sure
    #L value is too large weirdly
    
    #missing parameters that I need to convert into something somewhere
    #charge, mass? , dipole moment
    x0,y0,z0, vx0,vy0,vz0 = trsfrm.ctd2car(pitch, phase, Kinetic_energy, L_shell, latitude, longitude, mass, Re)
    #test = trsfrm.car2ctd(x0, y0, z0, vx0, vy0, vz0, mass, Re)
    #print(test)
    S0 = np.array([x0, y0, z0, vx0, vy0, vz0])
    S0 = S0/Re
    #print(S0)
    #non dimenionlize r. t is already nd
    #S0 = np.array([1,0.,0.,0.,.01,.01]) # temp override of initial conditions
    dT = 1/accuracy
    Bc = u0*M/(4*np.pi*Re**3)
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
    '''
    u0 = 1.25663706212e-6
    #B0 = M * u0 / (4*pi*Re^3) M = dipole moment
    M = 100
    B0 = M * u0 / (4*np.pi*Re**3)
    #tau = Mass / ( q * B0)
    tau = m / (e * B0)
    t = tau * T
    '''

    #particle_demo(pitch =j,_3dplots=False,fft = False, method = 'nd',_2dplots=(True))