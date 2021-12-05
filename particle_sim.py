
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
from matplotlib import pyplot as plt
from vtk_export import vtk_export
from joblib import Parallel
from joblib import delayed




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

#seperate plots and save


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
                 _2dplots=True,  # generate 2d plots or not
                 Save=False,  # when True saves large files with path to avoid
                 # repeated simulations
                 integrate=True,  # run integration again or load
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


    # choose integrator based of input
    if integrate is True:
        if method == 'rk45':
            xline, yline, zline, Vx, Vy, Vz, T = integrator.\
                rk45_nd(dT, tp, S0, qsign)
        elif method == 'euler':
            xline, yline, zline, Vx, Vy, Vz, T = integrator.\
                euler_cromer(dT, tp, S0, qsign)
        elif method == 'boris':
            xline, yline, zline, Vx, Vy, Vz, T,err_V = integrator.\
                boris(dT, sampling, S0, tp, qsign)
        else:
            print('invalid method')
        t = T*tc
        #print('integrator done at time ', datetime.now() - startTime)
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


def save(t, xline, yline, zline, Vx, Vy, Vz, L_shell, pitchangle,
         charge, mass, Kinetic_energy, method):
    if not os.path.exists('plots'):
        os.mkdir('plots')

    qm_rat = charge/mass
    kename = Kinetic_energy*1e-3

    filename = 'plots/qm_{:,}_Ke_{}MeV_pitch_{}d_L_{}Re_{}/'\
        .format(qm_rat, kename, pitchangle, L_shell, method)

    # title = '''pitch = {}\N{DEGREE SIGN} q/m = {:,}, L = {},
    #     Ke = {:.1e}eV, Method = {}''' \
    #     .format(pitchangle, qm_rat, L_shell, Kinetic_energy, method)
    # make directory for plots
    if not os.path.exists(filename):
        os.mkdir(filename)
    #tc = abs((mass/(charge*Bc)))
    #T = t/tc
    
    points = np.column_stack([xline[~np.isnan(xline)], yline[~np.isnan(yline)],
                              zline[~np.isnan(zline)]])
    out_filename = (filename + 'xyz.vtk')

    ftype = 'ASCII'
    connectivity = {'LINES': np.array([points.shape[0]])}

    vtk_export(out_filename, points,  # should points be lines?
               dataset='POLYDATA',
               connectivity=connectivity,
               title='Title',
               ftype=ftype,
               debug=False)

    arr_out = np.array((t, xline, yline, zline, Vx, Vy, Vz))
    np.savetxt(filename + 'trajectory.txt', arr_out)
    return


def plot(L_shell, pitchangle, charge, m, Kinetic_energy, method, err_V, units='s',
         t = None, xline = None, yline = None,zline = None,
         Vx = None, Vy = None, Vz = None):
    # valid choices for units are 's' 'Tc' 'min' 'days
    # x = None passes nothing if compute was not run, 
    # and pull from save if exists
    # but if a save exists, new compute get's priority for plot

    qm_rat = charge/m
    kename = Kinetic_energy*1e-3

    filename = 'plots/qm_{:,}_Ke_{}MeV_pitch_{}d_L_{}Re_{}/'\
        .format(qm_rat, kename, pitchangle, L_shell, method)

    title = '''pitch = {}\N{DEGREE SIGN} q/m = {:,}, L = {},
        Ke = {:.1e}eV, Method = {}''' \
        .format(pitchangle, qm_rat, L_shell, Kinetic_energy, method)
        
    # read in values if needed
    #check if file exists already?
    #how to check if other arrays exist?
    #if t == none then no array exists and need to read
    if t is None:
        if os.path.exists(filename):
            t, xline, yline, zline, Vx, Vy, Vz = np.loadtxt(
                filename+'trajectory.txt')
        elif not os.path.exists(filename):
            print('no such value was computed or saved')
            
    elif t is not None:
        #need to check if folder exists or not
        if not os.path.exists(filename):
            os.mkdir(filename)
    

    Ke = .5 * m * (np.power(Vx, 2) + np.power(Vy, 2) + np.power(Vz, 2))
    Ke0 = Ke[0]
    err_E = abs((Ke0-Ke)/Ke0)
    #  plot in seconds
    T = t
    if units == 's':
        timelabel = 'Time (seconds)'
        T_plot = t
    elif units == 'Tc':
        T_plot = T
        timelabel = 'Time (Tc)'
    elif units == 'min':
        T_plot = t/60
        timelabel = 'Time (minutes)'
    elif units == 'days':
        timelabel = 'Time (days)'
        T_plot = t/60/60/24
    else:
        print('invalid choice for units\nuse s, Tc, min, or days\n \
              defaulting to Tc')
        T_plot = T
        timelabel = 'Time (Tc)'
    # getting V parrallel and V perpendicular
    (Vparrallel, Vperpendicular, V_perp)\
        = trsfrm.VparVperp(xline, yline, zline, Vx, Vy, Vz)
    # for vtk export
    # xline[~np.isnan(xline) removes nan values, paraview can't work with nans

    # 2d Plotting: Energies, position, L-shell,x-y plot

    # title = 'pitch = {}, q/m ={}'.format(pitchangle, m)
    # energy plot
    

    size = [12.8, 9.6]
    plt.figure(1, figsize=size)
    plt.plot(T_plot,err_E) 
    plt.title(title)
    plt.xlabel(timelabel)
    plt.ylabel('Relative error in energy')
    plt.savefig(filename + 'energy.svg', format='svg')
    # plt.show()
    plt.close()
    plt.clf()
    # velocity plots
    # x,y,z

    # position
    plt.figure(4, figsize=size)
    plt.plot(T_plot, xline)
    # plt.plot(T_plot, yline)
    # plt.plot(T_plot, zline)
    plt.legend('x', loc='upper right')
    plt.title(title)
    plt.xlabel(timelabel)
    plt.ylabel('Distance (Re)')
    plt.savefig(filename + 'x.svg', format='svg')
    # plt.show()
    plt.close()
    plt.clf()

    plt.figure(5, figsize=size)
    # plt.plot(T_plot, xline)
    # plt.plot(T_plot, yline)
    plt.plot(T_plot, zline)
    plt.legend('z', loc='upper right')
    plt.title(title)
    plt.xlabel(timelabel)
    plt.ylabel('Distance (Re)')
    plt.savefig(filename + 'z.svg', format='svg')
    # plt.show()
    plt.close()
    plt.clf()

    plt.figure(7, figsize=size)
    # plt.plot(T_plot, xline)
    plt.plot(T_plot, yline)
    # plt.plot(T_plot, zline)
    plt.legend(['y'], loc='upper right')
    plt.title(title)
    plt.xlabel(timelabel)
    plt.ylabel('Distance (Re)')
    plt.savefig(filename + 'y.svg', format='svg')
    # plt.show()
    plt.close()
    plt.clf()
    # L-shell

    # xy plot
    plt.figure(6, figsize=size)
    plt.axes().set_aspect('equal', 'datalim')
    plt.plot(xline, yline)
    plt.title(title)
    plt.xlabel('X (Re)')
    plt.ylabel('Y (Re)')
    plt.savefig(filename + 'xy.svg', format='svg')
    # plt.show()
    plt.close()
    plt.clf()

def trajectory(pitch,m,Ke,q,T,acc,Lshell,sample,inte):
    t = 1e9*Ke**-1
    t, xline, yline, zline, Vx, Vy, Vz,err_V = particle_sim(
            pitchangle=pitch, mass=m, Kinetic_energy=Ke,
            charge=q, t=T, accuracy=acc, L_shell=Lshell,
            sampling=sample, method=inte, losscone=False)
    

    save(t, xline, yline, zline, Vx, Vy, Vz, Lshell, pitch, q,
          m, Ke, inte)

    plot(Lshell, pitch, q, m, Ke, inte,err_V, t=t, xline=xline, yline=yline,
          zline=zline, Vx=Vx, Vy=Vy, Vz=Vz)



def demo():
    m = [M_e, M_e, M_p, M_p, M_p, M_p, M_p, M_p]
    q = [C_e, C_e, -C_e, -C_e, -C_e, -C_e, -C_e, -C_e]
    Ke = [1e3, 1e3, 1e8, 1e8, 1e8, 1e8, 1e8, 1e8]
    T = [1e-5, 1, 1e-2, .2, 10, 10, 4, 4]
    acc = [1e2, 1e2, 1e2, 1e3, 1e2, 1e2, 1e4, 1e2]
    pitch = [90, 89, 90, 89, 90, 90, 5, 5]
    Lshell = [2.1, 1, 1.64, 1, 3, 3, 1, 1]
    sample = [20, .01, 10, 5, 1, 1, 1, 1]
    inte = ['boris', 'boris', 'boris', 'boris', 'boris', 'rk45',
            'boris', 'rk45']

    # takes my computer 30 mins to run all of these
    # wow that last one takes ages
    # why can't I parallize this? they share values so it comes out as nonsense
    
    Parallel(n_jobs=-2, prefer="threads")(
        delayed(trajectory)(pitch[i], m[i], Ke[i], q[i], T[i], 
        acc[i], Lshell[i], sample[i], inte[i]) for i in [0,2])
    
    
    # for i in range(0,4):
    #     trajectory(pitch[i], m[i], Ke[i], q[i], T[i],
    #                     acc[i], Lshell[i], sample[i], inte[i])
                       
        # len(m)):  # remeber to set this back to len(m)


def trajectory_generator():
    L = np.linspace(2, 10, 27)
    #pitch = np.linspace(90, 10, 8)
    pitch = 90
    m = [M_p, M_e]
    q = [-C_e, C_e]
    #K = np.logspace(8, 2, 7)
    K = 1e8
    T = 1
    acc= 1e1
# use joblib to speed up ranges
    # for a in range(0, len(m)):  # skip electron for now
    #     for b in range(0, len(K)):
    #         for c in range(0, len(pitch)):
    #             for d in range(0, len(L)):
    #                 # print(a+b+c+d)
    #                 tshell = 1e9*K[b]**-1
    #                 # print(tshell)
    #                 particle_sim(L_shell=L[d], pitchangle=pitch[c], mass=m[a],
    #                              charge=q[a],  Kinetic_energy=K[b], t=tshell)
    #                 #I think the longest ones will take around 8 days?
                    
    Parallel(n_jobs=-2, prefer='threads')(
            delayed(trajectory)(pitch[c], m[a], K[b], q[a], T, acc, 
             L[d],5, 'boris') for a in range(len(m)) for b in range(len(K))
            for c in range(len(pitch)) for d in range(len(L)))

demo()
#rajectory_generator()

