import numpy as np
from datetime import datetime

from joblib import Parallel
from joblib import delayed

import constants
import transformations as trsfrm
import integrator
from output import save
from output import plot
from particle_sim import particle_sim


def trajectory(pitch_angle, mass, Kinetic_energy, charge, accuracy, L_shell, sampling, method,t):

   

    #calculate estimated bouce and drift period
    R = L_shell # only valid when starting at equator
    pitch = np.radians(pitch_angle)
    # covert energy to joules
    
    Ke = Kinetic_energy * 1.602176565e-19
    phase,latitude,longitude = 0,0,0
    x0, y0, z0, vx0, vy0, vz0 = trsfrm.ctd2car(pitch,
                                               phase, Ke,
                                               L_shell, latitude, longitude,
                                               mass, constants.Re)
    v = np.linalg.norm([vx0,vy0,vz0])
    beta = v/constants.c
    #print(R,v,beta)
    tb = trsfrm.t_b(R,beta,pitch)
    #print('tb',tb)
    if mass == constants.M_p:
        Cd = 8.481
    elif mass == constants.M_e:
        Cd = 1.557e4
    else:
        print('drift calculation only for electrons and protons')
    td = trsfrm.t_d(R, beta, np.radians(pitch_angle),Cd)
    
   # print('bounce period is {:.2e} \ndrift period is {:.2e}'.format(tb,td))
    
    
    #print("Computing trajectory using " + method)
    T, xline, yline, zline, Vx, Vy, Vz = particle_sim(
        L_shell, pitch_angle, mass,charge, 1.1*td, Kinetic_energy, method, accuracy,
        sampling, losscone=False)

    # save(t, xline, yline, zline, Vx, Vy, Vz, Lshell, pitch, q,
    #       m, Ke, method)
    plot(L_shell, pitch_angle, charge, mass, Kinetic_energy, method, T,
          xline, yline, zline,
          Vx, Vy, Vz)

    
    
    

def trajectory_generator(par=True):
    L_shell = np.linspace(2, 10, 2)
    # pitch = np.linspace(90, 10, 8)
    pitch_angle = [90]
    mass = [constants.M_p,constants. M_e]
    charge = [-constants.C_e, constants.C_e]
    # K = np.logspace(8, 2, 7)
    Kinetic_energy = [1e8]
    
    #Kinetic_energy = Kinetic_energy * 1.602176565e-19
    accuracy = 1e2
    sampling = 36
    method = 'boris'
    
    if par is False:
        for a in range(0,1): # len(mass)):  # skip electron for now
            for b in range(0, len(Kinetic_energy)):
                for c in range(0, len(pitch_angle)):
                    for d in range(0, len(L_shell)):
                        trajectory(pitch_angle[c], mass[a], Kinetic_energy[b],
                                   charge[a], accuracy,L_shell[d], sampling, method,t = None)
                        


    # use joblib to speed up compute time
    # carful when running this, I think if memory gets cached, it will break
    # energy error plots
    if par is True:
        T = None
        Parallel(n_jobs=-2, prefer='threads')(
            delayed(trajectory)(pitch_angle[c], mass[a], Kinetic_energy[b],
                                 charge[a], accuracy, L_shell[d], sampling,
                                 method,T) 
                                for a in range(0, len(mass)) 
                                for b in range(len(Kinetic_energy))
                                for c in range(len(pitch_angle)) 
                                for d in range(len(L_shell))) 
        
trajectory_generator(par = True)
