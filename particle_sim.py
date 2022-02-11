
# This file perfroms all the transformations from the transformatons module
# computes the arrays using the Integrators module
# save to save plaintext and .vtk for paraview
# plot gernertes saved images to look at simualation parameters
import numpy as np
from datetime import datetime

from joblib import Parallel
from joblib import delayed

import transformations as trsfrm
import integrator
from output import save
from output import plot


startTime = datetime.now()
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


def particle_sim(L_shell=2,
                 pitch_angle=60,
                 mass=M_p,  # in Kg
                 charge=-C_e,  # in coulumbs
                 t=1e1,  # time in seconds
                 Kinetic_energy=1e8,
                 # in eV. Defaults to high energy to shorten drift period
                 
                 method='boris',# valid choices are 'boris','rk45'
                 # and 'euler__cromer'
                 accuracy=1e3,  # inverse time step in dimensionless form
                 sampling=30,   # points per gyro
                 # note accuracy*sampling cannot be greater than 1
                 losscone=True,# True to ditch atmoshperic particles
                 # False to keep them, like for demo()
                 latitude=0,  # all angles in degrees for input
                 longitude=0,
                 phase=0):
    

    
    #sampling modification
    sampling = sampling / (np.pi*2)
    
    sampling = sampling / L_shell**3
    accuracy = accuracy / L_shell**3
    
    # internally all angles are in radians
    latitude = np.radians(latitude)
    longitude = np.radians(longitude)
    phase = np.radians(phase)
    pitch = np.radians(pitch_angle)
    # covert energy to joules
    Kinetic_energy = Kinetic_energy * 1.602176565e-19

    dT = 1/accuracy
    # convert initial conditons to cartesian
    x0, y0, z0, vx0, vy0, vz0 = trsfrm.ctd2car(pitch,
                                               phase, Kinetic_energy,
                                               L_shell, latitude, longitude,
                                               mass, Re)
    S0 = np.array([x0, y0, z0, vx0, vy0, vz0])
    v = np.linalg.norm([vx0,vy0,vz0])
    beta = v/c
    #print(beta)
    gamma = (1-beta**2)**-.5
    mass = gamma*mass
    print(gamma)
    # convert into new dimensionless units
    # need constants first
    # re is at the top
    tc = abs((mass/(charge*Bc)))
    # in nd form can't have negative tc
    # qsign passes the sign of the charge
    # simply invert B if charge is negative.
    qsign = np.sign(charge)

    # convert to nd positions and velocities
    S0 = S0/Re
    S0[3:6] = S0[3:6] * tc
    # convert time into nd
    tp = t/(tc)
    # check loss cone angle and exit if angle is too high



    
    if losscone is True:
        B = np.linalg.norm(trsfrm.B(S0[0], S0[1], S0[2]))
        term = np.sqrt(B/Ba)
        if abs(term) > 1 or L_shell < 1:
            print('particle will hit the atmosphere, skipping')
            return

        alphac = np.arcsin(term)

        if pitch < alphac:
            print('particle will hit the atmosphere, skipping')
            return

    # choose integrator based of input
    if method == 'rk45':
        xline, yline, zline, Vx, Vy, Vz, T = integrator.\
            rk45_nd(dT, tp, S0, qsign)
    elif method == 'euler_cromer':
        xline, yline, zline, Vx, Vy, Vz, T = integrator.\
            euler_cromer(dT, tp, S0, qsign)
    elif method == 'boris':
        xline, yline, zline, Vx, Vy, Vz, T = integrator.\
            boris(dT, sampling, S0, tp, qsign)
    else:
        print('invalid method')
        return
    t = T*tc
        # print('integrator done at time ', datetime.now() - startTime)
    return t, xline, yline, zline, Vx, Vy, Vz


# information on particle perdiods can be found in period_info.pdf
# lists to generate various trajectory examples


def trajectory(pitch_angle, mass, Kinetic_energy, charge, t, accuracy, L_shell, sampling, method):

    if t is None:
        t = 1e9*Kinetic_energy

    t, xline, yline, zline, Vx, Vy, Vz = particle_sim(
        pitch_angle, mass, Kinetic_energy,
        charge, t, accuracy, L_shell,
        sampling, method, losscone=False)

    save(t, xline, yline, zline, Vx, Vy, Vz, L_shell, pitch_angle, charge,
         mass, Kinetic_energy, method)

    plot(L_shell, pitch_angle, charge, mass, Kinetic_energy, method, t, xline, yline,
         zline, Vx, Vy, Vz)

def trajectory_generator(par=True):
    L = np.linspace(2, 10, 7)
    # pitch = np.linspace(90, 10, 8)
    pitch = [90]
    m = [M_p, M_e]
    q = [-C_e, C_e]
    # K = np.logspace(8, 2, 7)
    K = [1e8]
    T = 1
    acc = 1e1

    if par is False:
        for a in range(0, len(m)):  # skip electron for now
            for b in range(0, len(K)):
                for c in range(0, len(pitch)):
                    for d in range(0, len(L)):
                        # print(a+b+c+d)
                        tshell = 1e9*K[b]**-1
                        # print(tshell)
                        particle_sim(L_shell=L[d], pitchangle=pitch[c],
                                     mass=m[a], charge=q[a],
                                     Kinetic_energy=K[b], t=tshell)

    # use joblib to speed up compute time
    # carful when running this, I think if memory gets chached, it will break
    # energy error plots
    if par is True:
        T = None
        Parallel(n_jobs=-2, prefer='threads')(
            delayed(trajectory)(pitch[c], m[a], K[b], q[a], T, acc,
                                L[d], 5, 'boris') for a in range(1, len(m))
            for b in range(len(K))
            for c in range(len(pitch)) for d in range(len(L)))


#demo()
# trajectory_generator()
