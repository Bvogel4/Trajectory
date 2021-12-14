import os
import numpy as np
from matplotlib import pyplot as plt

from vtk_export import vtk_export
import transformations as trsfrm

def save(t, xline, yline, zline, Vx, Vy, Vz, L_shell, pitchangle,
         charge, mass, Kinetic_energy, method):

    if not os.path.exists('output'):
        os.mkdir('output')

    qm_rat = charge/mass
    kename = Kinetic_energy*1e-3

    out_dir = 'output/qm_{:,}_Ke_{}MeV_pitch_{}d_L_{}Re_{}/'\
        .format(qm_rat, kename, pitchangle, L_shell, method)

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    # CSV
    out_filename = out_dir + 'trajectory.txt'

    arr_out = np.array((t, xline, yline, zline, Vx, Vy, Vz))

    # TODO: Output file should have 7 columns
    # Change velocities to R_E/s
    
    print("Saving " + out_filename)
    np.savetxt(out_filename, arr_out, header='t [s], x [R_E], y [R_E], z [R_E], vx [R_E/s], vy [R_E/s], vz [R_E/s]')
    print("Saved " + out_filename)


    # VTK
    out_filename = out_dir + 'xyz.vtk'

    # Boris method is passed nan-filled array and may not overwrite all fill
    # values.
    # TODO: Consider trimming output within Boris calculation.
    nans = np.isnan(xline)
    points = np.column_stack([xline[~nans], yline[~nans], zline[~nans]])

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


def plot(L_shell, pitchangle, charge, m, Kinetic_energy, method, err_V,
         units='s', t=None, xline=None, yline=None, zline=None, Vx=None,
         Vy=None, Vz=None):
    # valid choices for units are 's' 'Tc' 'min' 'days
    # x = None passes nothing if compute was not run,
    # and pull from save if exists
    # but if a save exists, new compute get's priority for plot

    if not os.path.exists('output'):
        os.mkdir('output')

    qm_rat = charge/m
    kename = Kinetic_energy*1e-3

    filename = 'output/qm_{:,}_Ke_{}MeV_pitch_{}d_L_{}Re_{}/'\
        .format(qm_rat, kename, pitchangle, L_shell, method)

    title = '''pitch = {}\N{DEGREE SIGN} q/m = {:,}, L = {},
        Ke = {:.1e}eV, Method = {}''' \
        .format(pitchangle, qm_rat, L_shell, Kinetic_energy, method)

    # read in values if needed
    # check if file exists already?
    # how to check if other arrays exist?
    # if t == none then no array exists and need to read
    if t is None:
        if os.path.exists(filename):
            t, xline, yline, zline, Vx, Vy, Vz = np.loadtxt(
                filename+'trajectory.txt')
        elif not os.path.exists(filename):
            print('no such value was computed or saved')

    elif t is not None:
        # need to check if folder exists or not
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
        print(
            'invalid choice for units use s, Tc, min, or days defaulting to s')
        timelabel = 'Time (seconds)'
        T_plot = t
    # getting V parrallel and V perpendicular
    (Vparrallel, Vperpendicular, V_perp)\
        = trsfrm.VparVperp(xline, yline, zline, Vx, Vy, Vz)
    # for vtk export
    # xline[~np.isnan(xline) removes nan values, paraview can't work with nans

    # 2d Plotting: Energies, position, L-shell,x-y plot

    # title = 'pitch = {}, q/m ={}'.format(pitchangle, m)
    # energy plot

    size = [12.8, 9.6]

    plt.figure(figsize=size)

    plt.plot(T_plot, err_E, label='error in energy')
    plt.title(title)
    plt.legend()
    plt.xlabel(timelabel)
    plt.ylabel('Relative error in energy')
    plt.savefig(filename + 'energy.svg', format='svg')
    plt.close()
    # plt.show()

    # velocity plots
    # x,y,z

    # position
    plt.figure(figsize=size)
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

    plt.figure(figsize=size)
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

    plt.figure(figsize=size)
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

    # L-shell

    # xy plot
    plt.figure(figsize=size)
    plt.axes().set_aspect('equal', 'datalim')
    plt.plot(xline, yline)
    plt.title(title)
    plt.xlabel('X (Re)')
    plt.ylabel('Y (Re)')
    plt.savefig(filename + 'xy.svg', format='svg')
    # plt.show()
    plt.close()
