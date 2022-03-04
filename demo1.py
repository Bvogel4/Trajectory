
import numpy as np
from joblib import Parallel
from joblib import delayed
from datetime import datetime


import transformations as trsfrm
from particle_sim import particle_sim
from output import save, plot


startTime = datetime.now()

#t, xline, yline, zline, Vx, Vy, Vz = None, None, None, None, None, None, None

#do a drift and bounce
def demo(mass,charge):

    #calculate estimated bouce and drift period
    R = L_shell # only valid when starting at equator
    pitch = np.radians(pitch_angle)
    # covert energy to joules
    Ke = Kinetic_energy * 1.602176565e-19
    phase,latitude,longitude = 0,0,0
    x0, y0, z0, vx0, vy0, vz0 = trsfrm.ctd2car(pitch,
                                               phase, Ke,
                                               L_shell, latitude, longitude,
                                               mass, Re)
    v = np.linalg.norm([vx0,vy0,vz0])
    beta = v/c
    #print(R,v,beta)
    tb = trsfrm.t_b(R,beta,pitch)
    #print('tb',tb)
    if mass == M_p:
        Cd = 8.481
    elif mass == M_e:
        Cd = 1.557e4
    else:
        print('drift calculation only for electrons and protons')
    td = trsfrm.t_d(R, beta, np.radians(pitch_angle),Cd)
    
    print('bounce period is {:.2e} \ndrift period is {:.2e}'.format(tb,td))
    
    
    #print("Computing trajectory using " + method)
    T, xline, yline, zline, Vx, Vy, Vz = particle_sim(
        L_shell, pitch_angle, mass,charge, 1.1*td, Kinetic_energy, method, accuracy,
        sampling, losscone=False)
    print('compute for {} done at time {}s'.format(method, datetime.now() - startTime))
    # save(t, xline, yline, zline, Vx, Vy, Vz, Lshell, pitch, q,
    #       m, Ke, method)
    plot(L_shell, pitch_angle, charge, mass, Kinetic_energy, method, T,
          xline, yline, zline,
          Vx, Vy, Vz)
    print('plot for {} done at time {}s'.format(method, datetime.now() - startTime))
    
    

    
    return

Re = 6.371e6
M_p = 1.6726219e-27  # kg
M_e = 9.10938356e-31  # kg
C_e = -1.60218e-19   # C
c = 3e8
mass = [M_p,M_e]
charge = [-C_e,C_e]
Kinetic_energy = 1e7    # eV
pitch_angle = 90  # degress
L_shell = 5
method = 'boris'
accuracy = 1e+2
sampling = 36

for a in range(len(mass)):
    demo(mass[a],charge[a])

#for method in methods:
    
    #lets do 1 drift and a few bounce for proton

# Parallel(n_jobs=1, prefer='threads')(delayed(demo)(L_shell[a],t[a])
#                                       for a  in range(len( L_shell)))
