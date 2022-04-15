import os
import numpy as np
from matplotlib import pyplot as plt

import constants
from vtk_export import vtk_export
import transformations as trsfrm


def save(parameters,t, x, y, z, vx, vy, vz):

    if not os.path.exists('output'):
        os.mkdir('output')

    species_name = parameters['species']
    
    if parameters['species']  != ('electron' or 'proton'):
        qm_rat = parameters['charge']/parameters['mass']
        species_name = 'qm_{:,}'.format(qm_rat)
    kename = parameters['Kinetic_energy']*1e-3

    out_dir = 'output/{}_Ke_{}MeV_pitch_{}d_L_{}Re_{}/'\
        .format(species_name, kename, parameters['pitch_angle'], 
                parameters['L_shell'], parameters['method'])

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    # CSV
    out_filename = out_dir + 'trajectory.csv'

    arr_out = np.column_stack((t, x, y, z, vx, vy, vz))

    # TODO: Output file should have 7 columns
    # Change velocities to R_E/s



    tc = abs((parameters['mass']/(parameters['charge']*constants.Bc)))

    vx = vx/tc
    vy = vy/tc
    vz = vz/tc

    #print("Saving " + out_filename)
    np.savetxt(out_filename, arr_out,
               header='t [s], x [R_E], y [R_E], z [R_E], vx [R_E/s], vy [R_E/s], vz [R_E/s]')
    #print("Saved " + out_filename)

    # VTK
    out_filename = out_dir + 'xyz.vtk'

    points = np.column_stack([x, y, z])

    ftype = 'ASCII'
    connectivity = {'LINES': np.array([points.shape[0]])}

    print("Saving " + out_filename)
    vtk_export(out_filename, points,  # should points be lines?
               dataset='POLYDATA',
               connectivity=connectivity,
               title='Title',
               ftype=ftype,
               debug=False)
    print("Saved " + out_filename)

    return


def plot(parameters,t=None, x=None, y=None, z=None, vx=None,vy=None, vz=None, units='s'):
    # valid choices for units are 's' 'Tc' 'min' 'days
    # x = None passes nothing if compute was not run,
    # and pull from save if exists
    # but if a save exists, new compute get's priority for plot

    if not os.path.exists('output'):
        os.mkdir('output')

    species_name = parameters['species']
    
    
    
    
    
    
    if not (parameters['species']   in ('electron' , 'proton')):
        qm_rat = parameters['charge']/parameters['mass']
        species_name = 'qm_{:,}'.format(qm_rat)
        
        
        
    kename = parameters['Kinetic_energy']*1e-3

    out_dir = 'output/{}_Ke_{}MeV_pitch_{}d_L_{}Re_{}/'\
        .format(species_name, kename, parameters['pitch_angle'], 
                parameters['L_shell'], parameters['method'])
        
        
    title = '''pitch = {}\N{DEGREE SIGN} q/m = {}, L = {},
        Ke = {:.1e}eV, Method = {}''' \
        .format(parameters['pitch_angle'], parameters['species'],
                parameters['L_shell'], parameters['Kinetic_energy'],
                parameters['method'])

    if t is None:
        if os.path.exists(out_dir):
            t, x, y, z, vx, vy, vz = np.loadtxt(
                out_dir+'trajectory.txt', unpack=True)
        elif not os.path.exists(out_dir):
            print('no such value was computed or saved')

    elif t is not None:
        # need to check if folder exists or not
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)

    #plot about 10 thousand points instead of billions
    length = len(t)
    if length > 10000:

        N = int(length/100000)

        t, x, y, z, vx, vy, vz, = t[::N], x[::N], y[::N], \
            z[::N], vx[::N], vy[::N], vz[::N]

    #error in energy
    Ke = .5 * parameters['mass'] * vx**2+vy**2+vz**2
    Ke0 = Ke[0]
    err_E = abs((Ke0-Ke)/Ke0)

    #  units for time
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
        print(
            'invalid choice for units use s, Tc, min, or days defaulting to s')
        timelabel = 'Time (seconds)'
        T_plot = t

    size = [12.8, 9.6]

    fig, ax = plt.subplots(figsize=size)

    ax.plot(T_plot, err_E, label='error in energy')
    ax.set_title(title)
    ax.legend()
    ax.set_xlabel(timelabel)
    ax.set_ylabel('Relative error in energy')
    fig.savefig(out_dir + 'energy.svg', format='svg')
    plt.close()
    # ax.show()

    # velocity plots
    # x,y,z

    # position
    fig, ax = plt.subplots(figsize=size)
    ax.plot(T_plot, x)
    # ax.plot(T_plot, y)
    # ax.plot(T_plot, z)
    ax.legend('x', loc='upper right')
    ax.set_title(title)
    ax.set_xlabel(timelabel)
    ax.set_ylabel('Distance (Re)')
    fig.savefig(out_dir + 'x.svg', format='svg')
    # ax.show()
    plt.close()

    fig, ax = plt.subplots(figsize=size)
    # ax.plot(T_plot, x)
    # ax.plot(T_plot, y)
    ax.plot(T_plot, z)
    ax.legend('z', loc='upper right')
    ax.set_title(title)
    ax.set_xlabel(timelabel)
    ax.set_ylabel('Distance (Re)')
    fig.savefig(out_dir + 'z.svg', format='svg')
    # ax.show()
    plt.close()

    fig, ax = plt.subplots(figsize=size)
    # ax.plot(T_plot, x)
    ax.plot(T_plot, y)
    # ax.plot(T_plot, z)
    ax.legend(['y'], loc='upper right')
    ax.set_title(title)
    ax.set_xlabel(timelabel)
    ax.set_ylabel('Distance (Re)')
    fig.savefig(out_dir + 'y.svg', format='svg')
    # ax.show()
    plt.close()

    # L-shell

    # xy plot
    fig, ax = plt.subplots(figsize=size)
    ax.set_aspect('equal', 'datalim')
    ax.plot(x, y)
    ax.set_title(title)
    ax.set_xlabel('X (Re)')
    ax.set_ylabel('Y (Re)')
    fig.savefig(out_dir + 'xy.svg', format='svg')
    # ax.show()
    plt.close()
