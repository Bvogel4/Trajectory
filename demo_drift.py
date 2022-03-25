import numpy as np
from joblib import Parallel
from joblib import delayed
from datetime import datetime

import constants
import transformations as trsfrm
from particle_sim import particle_sim
from output import save, plot


startTime = datetime.now()

#t, xline, yline, zline, Vx, Vy, Vz = None, None, None, None, None, None, None

#do a drift and bounce


def demo(mass, charge):

    #calculate estimated bouce and drift period
    R = L_shell  # only valid when starting at equator
    pitch = np.radians(pitch_angle)
    # covert energy to joules
    Ke = Kinetic_energy * 1.602176565e-19
    phase, latitude, longitude = 0, 0, 0
    x0, y0, z0, vx0, vy0, vz0 = trsfrm.ctd2car(pitch,
                                               phase, Ke,
                                               L_shell, latitude, longitude,
                                               mass, constants.Re)
    v = np.linalg.norm([vx0, vy0, vz0])
    beta = v/constants.c
    #print(R,v,beta)
    tb = trsfrm.t_b(R, beta, pitch)
    #print('tb',tb)
    if mass == constants.M_p:
        Cd = 8.481
    elif mass == constants.M_e:
        Cd = 1.557e4
    else:
        print('drift calculation only for electrons and protons')
    td = trsfrm.t_d(R, beta, np.radians(pitch_angle), Cd)

    print('bounce period is {:.2e} \ndrift period is {:.2e}'.format(tb, td))

    #print("Computing trajectory using " + method)
    T, xline, yline, zline, Vx, Vy, Vz = particle_sim(
        L_shell, pitch_angle, mass, charge, 1*td, Kinetic_energy, method, accuracy,
        sampling, losscone=False)
    print('compute for {} done at time {}s'.format(
        method, datetime.now() - startTime))
    # save(t, xline, yline, zline, Vx, Vy, Vz, Lshell, pitch, q,
    #       m, Ke, method)
    plot(L_shell, pitch_angle, charge, mass, Kinetic_energy, method, T,
         xline, yline, zline,
         Vx, Vy, Vz)
    print('plot for {} done at time {}s'.format(
        method, datetime.now() - startTime))

    return


mass = [constants.M_p, constants.M_e]
charge = [-constants.C_e, constants.C_e]
Kinetic_energy = 1e7    # eV
pitch_angle = 90  # degress
L_shell = 5
method = 'boris'
accuracy = 1e+2
sampling = 36

for a in range(len(mass)):
    demo(mass[a], charge[a])

#for method in methods:

    #lets do 1 drift and a few bounce for proton

# Parallel(n_jobs=1, prefer='threads')(delayed(demo)(L_shell[a],t[a])
#                                       for a  in range(len( L_shell)))
