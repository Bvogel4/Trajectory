# This file perfroms all the transformations from the transformatons module
# computes the arrays using the Integrators module
# save to save plaintext and .vtk for paraview
# plot gernertes saved images to look at simualation parameters
import numpy as np
from datetime import datetime

import transformations as trsfrm
import integrator
import constants

startTime = datetime.now()


def particle_sim(L_shell=2,
                 pitch_angle=60,
                 mass=constants.M_p,  # in Kg
                 charge=-constants.C_e,  # in coulumbs
                 t=1e1,  # time in seconds
                 Kinetic_energy=1e8,
                 # in eV. Defaults to high energy to shorten drift period

                 method='boris',  # valid choices are 'boris','rk45'
                 # and 'euler__cromer'
                 accuracy=1e3,  # inverse time step in dimensionless form
                 sampling=30,   # points per gyro
                 # note accuracy*sampling cannot be greater than 1
                 losscone=True,  # True to ditch atmoshperic particles
                 # False to keep them, like for demo()
                 latitude=0,  # all angles in degrees for input
                 longitude=0,
                 phase=0):

    #print(L_shell,pitch_angle,mass,charge,t,Kinetic_energy,method,accuracy,
    #sampling,losscone,latitude,longitude,phase)

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
                                               mass, constants.Re)
    S0 = np.array([x0, y0, z0, vx0, vy0, vz0])
    v = np.linalg.norm([vx0, vy0, vz0])
    beta = v/constants.c

    gamma = (1-beta**2)**-.5
    mass = gamma*mass

    # convert into new dimensionless units
    # need constants first
    # re is at the top
    tc = abs((mass/(charge*constants.Bc)))
    # in nd form can't have negative tc
    # qsign passes the sign of the charge
    # simply invert B if charge is negative.
    qsign = np.sign(charge)

    # convert to nd positions and velocities
    S0 = S0/constants.Re
    S0[3:6] = S0[3:6] * tc
    # convert time into nd
    tp = t/(tc)

    # check loss cone angle and exit if angle is too high
    if losscone is True:
        B = np.linalg.norm(trsfrm.B(S0[0], S0[1], S0[2]))
        term = np.sqrt(B/constants.Ba)
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

        # print('integrator done at time ', datetime.now() - startTime)
    # convert time into seconds
    t = T*tc
    #convert velocites into Re/s
    Vx = Vx / tc
    Vy = Vy / tc
    Vz = Vz / tc

    return t, xline, yline, zline, Vx, Vy, Vz

# information on particle perdiods can be found in period_info.pdf
