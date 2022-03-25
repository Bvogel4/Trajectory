
import numpy as np
from datetime import datetime
from datetime import timedelta

import constants
import transformations as trsfrm
from particle_sim import particle_sim
from output import plot


#t, xline, yline, zline, Vx, Vy, Vz = None, None, None, None, None, None, None


def demo_method(method, pitch_angle, accuracy):
    startTime = datetime.now()
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
    #print(beta)
    #print(R,v,beta)
    #tb = trsfrm.t_b(R, beta, np.radians(pitch_angle))
    #print('tb',tb)
    if mass == constants.M_p:
        Cd = 8.481
        td = trsfrm.t_d(R, beta, np.radians(pitch_angle), Cd)
    elif mass == constants.M_e:
        Cd = 1.557e4
        td = trsfrm.t_d(R, beta, np.radians(pitch_angle), Cd)
    else:
        print('drift calculation only for electrons and protons')

    #print('bounce period is {:.2e} \ndrift period is {:.2e}'.format(tb,td))

    #print("Computing trajectory using " + method)
    T, xline, yline, zline, Vx, Vy, Vz = particle_sim(
        L_shell, pitch_angle, mass, charge, td, Kinetic_energy, method,
        accuracy, sampling, losscone=False)
    compute_time = timedelta.total_seconds(datetime.now() - startTime)
    print('compute for {} done at time {}s'.format(
        method, compute_time))
    # save(t, xline, yline, zline, Vx, Vy, Vz, Lshell, pitch, q,
    #       m, Ke, method)
    plot(L_shell, pitch_angle, charge, mass, Kinetic_energy, method, T,
         xline, yline, zline,
         Vx, Vy, Vz)
    # print('plot for {} done at time {}s'.format(method,
    #                                             datetime.now() - startTime))

    Ke0 = (np.power(Vx[0], 2) + np.power(Vy[0], 2) + np.power(Vz[0], 2))
    Kef = (np.power(Vx[-1], 2) + np.power(Vy[-1], 2) + np.power(Vz[-1], 2))
    err_E = abs((Kef-Ke0)/Ke0)
    compute_efficiency = err_E/compute_time
    print(compute_efficiency)

    return compute_efficiency


mass = constants.M_p
charge = -constants.C_e
Kinetic_energy = 1e8    # eV

L_shell = 1

accuracy = [1e2, 1e4, 1e2, 1e2, 1e4, 1e2]
sampling = 36

method = ['boris', 'boris', 'rk45', 'boris', 'boris', 'rk45']
pitch_angle = [90, 90, 90, 10, 10, 10]
compute_efficiency = np.zeros(6)*np.nan 

for a in [3,4]:#range(len(pitch_angle)):
    compute_efficiency[a] = demo_method(method[a], pitch_angle[a], accuracy[a])

# Parallel(n_jobs=-2, prefer='threads')(delayed(demo_method)(methods[a],pitch[a])
#                                       for a in range(len(pitch)))
