# -*- coding: utf-8 -*-
"""
Created on Sun Aug 22 17:21:12 2021

@author: blake
"""
import numpy as np
import Transformations as trsfrm
import Integrator as Integrator
import plots
Re = 1 #6.37e6 # radius of earth in m
M = 100 #magnetic moment M = 8x10^22 \A m^2 for earth
u0 = u0 = 1.25663706212e-6
#for j in [90,60,45,30,0]:
for pitchangle in [90]: # generates plots with different intial pitch angles
    filename = 'pitch{0}'.format(pitchangle)
    L_shell = 1 
    latitude = 0 #all angles in degrees
    longitude = 0 
    phase = 0 
    Kinetic_energy = 1  # in eV/kg
    mass = 1e-18  #in Kg
    charge = 5e-18  # in coulumbs
    tp = 1 #time in tc ~ 50 error gets high for rk45, boris is much more accurate
    accuracy = 1e4
    method = 'rk45'
    _2dplots = True
    
    
    #convert initial conditons to cartesian
    #is ctd2car working right? I don't think it is
    #!!!! need to fix, but need to fix integrator first
    #think Re is a missing paramter
    x0,y0,z0, vx0,vy0,vz0 = trsfrm.ctd2car(pitchangle, phase, Kinetic_energy, L_shell, latitude, longitude, mass,Re)
    test = trsfrm.car2ctd(x0, y0, z0, vx0, vy0, vz0, mass, Re)
    print(test)
    S0 = np.array([x0, y0, z0, vx0, vy0, vz0])
    
    #print(S0)
    #non dimenionlize r. t is already nd
    if method == 'rk45':
        Bc = u0*M/(4*np.pi*Re**3)
        #tc = mass/(charge*Bc)
        
        dT = tp/accuracy
        #print(Bc,dT)
        #nd initial conditions
        S0[0:3] = S0[0:3]/Re
        S0 = [1.,0.,0.,1.,0.,0.]
        xline,yline,zline,Vx,Vy,Vz,T = Integrator.rk45_nd(dT, tp, S0)
        #convert back into dimensional coords?
        
        #print(T)
        #print(xline,yline)
        
        if _2dplots == True:
            tc = mass/(charge*Bc)
            t = T*tc
            plots.plotter(xline,yline,zline,Vx,Vy,Vz,filename, T,method,pitchangle,mass,Re)
'''
    elif method == 'boris':
        casdf = 2
    
    else:
        print('invalid method')
'''
    
    
    
    #particle_demo(pitch =j,_3dplots=False,fft = False, method = 'nd',_2dplots=(True))