import os
import numpy as np
from matplotlib import pyplot as plt

from vtk_export import vtk_export
import transformations as trsfrm

def save(t, xline, yline, zline, Vx, Vy, Vz, L_shell, pitch_angle,
         charge, mass, Kinetic_energy, method):

    if not os.path.exists('output'):
        os.mkdir('output')

    qm_rat = charge/mass
    kename = Kinetic_energy*1e-3

    out_dir = 'output/qm_{:,}_Ke_{}MeV_pitch_{}d_L_{}Re_{}/'\
        .format(qm_rat, kename, pitch_angle, L_shell, method)

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    # CSV
    out_filename = out_dir + 'trajectory.csv'

    arr_out = np.column_stack((t, xline, yline, zline, Vx, Vy, Vz))

    # TODO: Output file should have 7 columns
    # Change velocities to R_E/s
    
    
    Re = 6.371e6
    # Re = 1 # radius of earth in m
    M = 8.22e22
    # magnetic moment M = 8.22x10^22 \A m^2 for earth
    u0 = 1.25663706212e-6
    # critical B value
    Bc = u0*M/(4*np.pi*Re**3)
    tc = abs((mass/(charge*Bc)))
    
    Vx = Vx/tc
    Vy = Vy/tc
    Vz = Vz/tc
    
    #print("Saving " + out_filename)
    np.savetxt(out_filename, arr_out, header=
        't [s], x [R_E], y [R_E], z [R_E], vx [R_E/s], vy [R_E/s], vz [R_E/s]')
    #print("Saved " + out_filename)


    # VTK
    out_filename = out_dir + 'xyz.vtk'

    points = np.column_stack([xline, yline, zline])

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


def plot(L_shell, pitch_angle, charge, mass, Kinetic_energy, method,
          t=None, xline=None, yline=None, zline=None, Vx=None,
         Vy=None, Vz=None, units='s'):
    # valid choices for units are 's' 'Tc' 'min' 'days
    # x = None passes nothing if compute was not run,
    # and pull from save if exists
    # but if a save exists, new compute get's priority for plot

    if not os.path.exists('output'):
        os.mkdir('output')

    qm_rat = charge/mass
    kename = Kinetic_energy*1e-3

    filename = 'output/qm_{:,}_Ke_{}MeV_pitch_{}d_L_{}Re_{}/'\
        .format(qm_rat, kename, pitch_angle, L_shell, method)

    title = '''pitch = {}\N{DEGREE SIGN} q/m = {:,}, L = {},
        Ke = {:.1e}eV, Method = {}''' \
        .format(pitch_angle, qm_rat, L_shell, Kinetic_energy, method)

    # read in values if needed
    # check if file exists already?
    # how to check if other arrays exist?
    # if t == none then no array exists and need to read
    if t is None:
        if os.path.exists(filename):
            t, xline, yline, zline, Vx, Vy, Vz = np.loadtxt(
                filename+'trajectory.txt',unpack = True)
        elif not os.path.exists(filename):
            print('no such value was computed or saved')

    elif t is not None:
        # need to check if folder exists or not
        if not os.path.exists(filename):
            os.mkdir(filename)
            
    #plot about a thousands points instead of billions
    length = len(t)
    if length > 5000:
        
        N = int(length/5000)

        t, xline, yline, zline, Vx, Vy, Vz, = t[::N], xline[::N], yline[::N], \
            zline[::N], Vx[::N], Vy[::N], Vz[::N]

    #error in energy
    Ke = .5 * mass * (np.power(Vx, 2) + np.power(Vy, 2) + np.power(Vz, 2))
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
    

    


    fig,ax = plt.subplots(figsize=size)

    ax.plot(T_plot, err_E, label='error in energy')
    ax.set_title(title)
    ax.legend()
    ax.set_xlabel(timelabel)
    ax.set_ylabel('Relative error in energy')
    fig.savefig(filename + 'energy.svg', format='svg')
    plt.close()
    # ax.show()

    # velocity plots
    # x,y,z

    # position
    fig,ax = plt.subplots(figsize=size)
    ax.plot(T_plot, xline)
    # ax.plot(T_plot, yline)
    # ax.plot(T_plot, zline)
    ax.legend('x', loc='upper right')
    ax.set_title(title)
    ax.set_xlabel(timelabel)
    ax.set_ylabel('Distance (Re)')
    fig.savefig(filename + 'x.svg', format='svg')
    # ax.show()
    plt.close()

    fig,ax = plt.subplots(figsize=size)
    # ax.plot(T_plot, xline)
    # ax.plot(T_plot, yline)
    ax.plot(T_plot, zline)
    ax.legend('z', loc='upper right')
    ax.set_title(title)
    ax.set_xlabel(timelabel)
    ax.set_ylabel('Distance (Re)')
    fig.savefig(filename + 'z.svg', format='svg')
    # ax.show()
    plt.close()

    fig,ax = plt.subplots(figsize=size)
    # ax.plot(T_plot, xline)
    ax.plot(T_plot, yline)
    # ax.plot(T_plot, zline)
    ax.legend(['y'], loc='upper right')
    ax.set_title(title)
    ax.set_xlabel(timelabel)
    ax.set_ylabel('Distance (Re)')
    fig.savefig(filename + 'y.svg', format='svg')
    # ax.show()
    plt.close()

    # L-shell

    # xy plot
    fig,ax = plt.subplots(figsize=size)
    ax.set_aspect('equal', 'datalim')
    ax.plot(xline, yline)
    ax.set_title(title)
    ax.set_xlabel('X (Re)')
    ax.set_ylabel('Y (Re)')
    fig.savefig(filename + 'xy.svg', format='svg')
    # ax.show()
    plt.close()
