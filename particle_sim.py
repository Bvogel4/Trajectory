
# This file perfroms all the transformations from the transformatons module
# computes the arrays using the Integrators module
# save to save plaintext and .vtk for paraview
# plot gernertes saved images to look at simualation parameters
import math
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
# critical B value
Bc = u0*M/(4*np.pi*Re**3)

M_p = 1.6726219e-27  # kg
M_e = 9.10938356e-31
C_e = -1.60218e-19  # C
# magnetic field at the top of the atmopshere in nd terms
za = (Re + 1e5)/Re
Ba = np.linalg.norm(trsfrm.B(0, 0, za))


def particle_sim(L_shell=2,
                 pitchangle=60,
                 latitude=0,  # all angles in degrees for input
                 longitude=0,
                 phase=0,
                 Kinetic_energy=1e8,
                 # in eV. Defaults to high energy to shorten drift period
                 # default model = alpha particle
                 mass=M_p,  # in Kg
                 charge=-C_e,  # in coulumbs
                 t=1e1,  # time in seconds
                 # 1 tp corresponds to the time to go 1 radian
                 # around the gyro radius at L = 1
                 accuracy=1e3,  # inverse time step in dimensionless form
                 sampling=5,  # points per tc
                 # note accuracy*sampling cannot be greater than 1
                 method='boris',  # valid choices are 'boris','rk45'
                                  # and 'euler'
                 # previously generated trajectories.
                 losscone=True):  # True to ditch atmoshperic particles
    # False to keep them, like for demo()

    # internally all angles are in radians
    latitude = math.radians(latitude)
    longitude = math.radians(longitude)
    phase = math.radians(phase)
    pitch = math.radians(pitchangle)
    # covert energy to joules
    Kinetic_energy = Kinetic_energy * 1.602176565e-19

    dT = 1/accuracy
    # convert initial conditons to cartesian
    x0, y0, z0, vx0, vy0, vz0 = trsfrm.ctd2car(pitch,
                                               phase, Kinetic_energy,
                                               L_shell, latitude, longitude,
                                               mass, Re)
    S0 = np.array([x0, y0, z0, vx0, vy0, vz0])
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

    err_V = None
    # choose integrator based of input
    if method == 'rk45':
        xline, yline, zline, Vx, Vy, Vz, T = integrator.\
            rk45_nd(dT, tp, S0, qsign)
    elif method == 'euler':
        xline, yline, zline, Vx, Vy, Vz, T = integrator.\
            euler_cromer(dT, tp, S0, qsign)
    elif method == 'boris':
        xline, yline, zline, Vx, Vy, Vz, T, err_V = integrator.\
            boris(dT, sampling, S0, tp, qsign)
    else:
        print('invalid method')
    t = T*tc

    if type(err_V) != np.ndarray:

        if err_V is None:
            V = np.linalg.norm(np.array((Vx, Vy, Vz)), axis=0)
            err_V = abs(V[0] - V)/V[0]
        # print('integrator done at time ', datetime.now() - startTime)
        '''
        if Save is True:
            # save to text
            # v = v*tc to get in re/s
            arr_out = np.column_stack((t, xline, yline, zline, Vx, Vy, Vz))
            # add comment for savetext
            np.savetxt(filename + 'trajectory.txt', arr_out)
            # colums are time seconds, x,y,z in Vx, in re/s
    # load data
    if integrate is False:
        trajectory = np.loadtxt(filename+'trajectory.txt')
        t, xline, yline, zline, Vx, Vy, Vz = trajectory.T

    if _2dplots is True:
        # Bc = u0*M/(4*np.pi*Re**3)
        # tc = abs(mass/(charge*Bc))
        # t = T*tc
        plots.plotter(xline, yline, zline, Vx, Vy, Vz, filename,
                      title, T, t, method, pitchangle, mass, Re, 's')
        print('plotter done at time ', datetime.now() - startTime)
        '''
    return t, xline, yline, zline, Vx, Vy, Vz, err_V


# information on particle perdiods can be found in period_info.pdf
# lists to generate various trajectory examples


def trajectory(pitch, m, Ke, q, T_final, acc, Lshell, sample, inte):

    if T_final is None:
        T_final = 1e9*Ke**-1

    t, xline, yline, zline, Vx, Vy, Vz, err_V = particle_sim(
        pitchangle=pitch, mass=m, Kinetic_energy=Ke,
        charge=q, t=T_final, accuracy=acc, L_shell=Lshell,
        sampling=sample, method=inte, losscone=False)

    save(t, xline, yline, zline, Vx, Vy, Vz, Lshell, pitch, q,
         m, Ke, inte)

    plot(Lshell, pitch, q, m, Ke, inte, err_V, t=t, xline=xline, yline=yline,
         zline=zline, Vx=Vx, Vy=Vy, Vz=Vz)

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
