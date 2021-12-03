
# This file perfroms all the transformations from the transformatons module
# computes the arrays using the Integrators module
# plots using the plots module
import math
import os
import numpy as np
import transformations as trsfrm
import integrator
import plots
from datetime import datetime

if not os.path.exists('plots'):
    os.mkdir('plots')

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
                 accuracy=1e4,  # inverse time step in dimensionless form
                 sampling=1,  # points per tc
                 # note accuracy*sampling cannot be greater than 1
                 method='boris',  # valid choices are 'boris','rk45'
                                  # and 'euler'
                 _2dplots=True,  # generate 2d plots or not
                 Save=False,  # when True saves large files with path to avoid
                 # repeated simulations
                 integrate=True,  # run integration again or load
                 # previously generated trajectories.
                 losscone=True):  # True to ditch atmoshperic particles
    # False to keep them, like for demo()

    mname = (mass/M_p)
    qname = abs(int((charge/C_e)))
    qm_rat = charge/mass
    kename = Kinetic_energy*1e-3
    filename = 'plots/qm_{:,}_Ke_{}MeV_pitch_{}d_L_{}Re_{}/'\
        .format(qm_rat, kename, pitchangle, L_shell, method)

    title = '''pitch = {}\N{DEGREE SIGN} ,m = {:.2}amu,q = {}e,L = {},
        Ke = {:1e}eV, Method = {}''' \
        .format(pitchangle, mname, qname, L_shell, Kinetic_energy, method)
    # make directory for plots
    if not os.path.exists(filename):
        os.mkdir(filename)

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
        if abs(term) > 1:
            print('particle will hit the atmosphere, skipping')
            return

    alphac = np.arcsin(term)
    # print(pitch)
    # print(alphac)
    # there is a missed range when L is less than Re
    if pitch < alphac:
        print('particle will hit the atmosphere, skipping')
        return
    # sampling = 10/t

    # choose integrator based of input
    if integrate is True:
        if method == 'rk45':
            xline, yline, zline, Vx, Vy, Vz, T = integrator.\
                rk45_nd(dT, tp, S0, qsign)
        elif method == 'euler':
            xline, yline, zline, Vx, Vy, Vz, T = integrator.\
                euler_cromer(dT, tp, S0, qsign)
        elif method == 'boris':
            xline, yline, zline, Vx, Vy, Vz, T = integrator.\
                boris(dT, sampling, S0, tp, qsign)
        else:
            print('invalid method')

        print('integrator done at time ', datetime.now() - startTime)
        if Save is True:
            # save to text
            t = T*tc
            # v = v*tc to get in re/s
            arr_out = np.column_stack((t, xline, yline, zline, Vx, Vy, Vz))
            # add comment for savetext
            np.savetxt(filename + 'trajectory.txt', arr_out)
            # colums are time seconds, x,y,z in Vx, in re/s
    print('integrator done at time ', datetime.now() - startTime)
    # load data
    if integrate is False:
        trajectory = np.loadtxt(filename+'trajectory.txt')
        xline, yline, zline, Vx, Vy, Vz, T = trajectory.T

    if _2dplots is True:
        # Bc = u0*M/(4*np.pi*Re**3)
        # tc = abs(mass/(charge*Bc))
        t = T*tc
        plots.plotter(xline, yline, zline, Vx, Vy, Vz, filename,
                      title, T, t, method, pitchangle, mass, Re, 's')
        print('plotter done at time ', datetime.now() - startTime)


# information on particle perdiods can be found in period_info.pdf
# lists to generate various trajectory examples


def demo():
    m = [M_e, M_e, M_p, M_p, M_p, M_p, M_p, M_p]
    q = [C_e, C_e, -C_e, -C_e, -C_e, -C_e, -C_e, -C_e]
    Ke = [1e3, 1e3, 1e8, 1e8, 1e8, 1e8, 1e8, 1e8]
    T = [1e-5, 1, 1e-2, .2, 10, 10, 10, 10]
    acc = [1e2, 1e2, 1e2, 1e3, 1e2, 1e2, 1e4, 1e2]
    pitch = [90, 89, 90, 89, 90, 90, 5, 5]
    Lshell = [2.1, 1, 1.64, 1, 3, 3, 1, 1]
    sample = [10, .01, 10, 5, 1, 1, 1, 1]
    inte = ['boris', 'boris', 'boris', 'boris', 'boris', 'rk45',
            'boris', 'rk45']
    # loss cone angle

    # compute what sampling should be to have reasoonable file sizes?
    # shows 5 cases showing different trajectores, particles, and integrators
    # takes my computer ~ mins to run all of these

    for i in range(0, len(m)):  # remeber to se this back to 0
        particle_sim(pitchangle=pitch[i], mass=m[i], Kinetic_energy=Ke[i],
                     charge=q[i], t=T[i], accuracy=acc[i], L_shell=Lshell[i],
                     sampling=sample[i], method=inte[i], losscone=False)


demo()
# particle_sim(pitchangle=90, L_shell=2)


def trajectory_generator():  # this will probably take forever to run
    L = np.linspace(2, 10, 7)
    m = [M_p, M_e]
    pitch = np.linspace(90, 10, 8)
    K = np.logspace(2, 8, 7)

    for a in range(0, 1):  # len(m)):  # skip electron for now
        for b in range(0, 1):  # len(K)):
            for c in range(0, 1):  # len(pitch)):
                for d in range(1, len(L)):
                    tshell = 1e9*K[b]**-1
                    particle_sim(
                        L_shell=L[d], pitchangle=pitch[c], mass=m[a],
                        Kinetic_energy=K[b], t=tshell)


# trajectory_generator()
