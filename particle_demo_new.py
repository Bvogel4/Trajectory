
# This file perfroms all the transformations from the transformatons module
# computes the arrays using the Integrators module
# plots using the plots module
from datetime import datetime

# from numba import njit, prange
import numpy as np
import transformations as trsfrm
import integrator
import plots
import math
import os

startTime = datetime.now()
# constants
Re = 6.371e6
# Re = 1 # radius of earth in m
M = 8.22e22
# magnetic moment M = 8.22x10^22 \A m^2 for earth
u0 = 1.25663706212e-6

M_p = 1.6726219e-27  # kg
M_e = 9.10938356e-31
C_e = -1.60218e-19  # C


if not os.path.exists('particle_plots'):
    os.mkdir('particle_plots')


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
                 t=1e3,  # time in seconds
                 # 1 tp corresponds to the time to go 1 radian
                 # around the gyro radius at L = 1
                 accuracy=1e3,  # inverse time step in dimensionless form
                 sampling=5,  # points per tc
                 # note accuracy*sampling cannot be greater than 1
                 method='boris',  # valid choices are 'boris','rk45'
                                  # and 'euler'
                 _2dplots=True,  # generate 2d plots or not
                 integrate=True):  # run integration again or load
                                   # previously generated trajectories.

    mname = (mass/M_p)
    qname = abs(int((charge/C_e)))

    filename = 'particle_plots/pitch_{}_m_{:.2}amu_q_{:}e_L_{}/'\
        .format(pitchangle, mname, qname, L_shell)

    title = 'pitch = {}\N{DEGREE SIGN} ,m = {:.2}amu,q = {}e, L = {}, Ke = {:1e}eV' \
        .format(pitchangle, mname, qname, L_shell, Kinetic_energy)
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
    Bc = u0*M/(4*np.pi*Re**3)
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
        # save to text
        arr_out = np.column_stack((xline, yline, zline, Vx, Vy, Vz, T))
        np.savetxt(filename + 'trajectory.txt', arr_out)
        print('integrator done at time ', datetime.now() - startTime)
        # save other useful info for xxxx?

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


# iterate over 2 protons, and 1 alpha particle
# electrons may need to be done alone, very memory intensive
# information on particle perdiods can be found in period_info.pdf
# lists to generate various trajectory examples
m = [M_e, M_e, M_p, M_p, M_p]
q = [C_e, C_e, -C_e, -C_e, -C_e]
Ke = [1e3, 1e3, 1e8, 1e8, 1e8]
T = [1e-5, 1, 1e-2, .2, 10]
acc = [1e2, 1e2, 1e2, 1e3, 1e2]
pitch = [90, 89, 90, 89, 90]
Lshell = [2.1, 1, 1.64, 1, 3]
sample = [10, .1, 10, 5, 1]
inte = ['boris', 'boris', 'boris', 'boris', 'boris']
# name = ['equatorial_proton_L2','Proton_L2_pitch.60','Equatorial_electron_L8',
#        'equatroial_alpha_L2','alpha_L2_pitch.60','alpha_L2_pitch.30']
# particle_sim()

# shows 5 cases showing different trajectores, particles, and integrators
# takes my computer ~8 mins to run all of these


for i in range(len(m)):
    particle_sim(pitchangle=pitch[i], mass=m[i], Kinetic_energy=Ke[i],
                 charge=q[i], t=T[i], accuracy=acc[i],
                 L_shell=Lshell[i], sampling=sample[i], method=inte[i])
