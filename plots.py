
import numpy as np
from matplotlib import pyplot as plt
from vtk_export import vtk_export
import transformations as trsfrm


def plotter(xline, yline, zline, Vx, Vy, Vz, filename, title,
            T, t,  # nd and seconds time
            integrator, pitchangle, m, Re, units):
    # valid choices for units are 's' 'Tc' 'min' 'days
    Ke = .5 * m * (np.power(Vx, 2) + np.power(Vy, 2) + np.power(Vz, 2))
    Ke0 = Ke[0]
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
    # 2d Plotting: Energies, position, L-shell,x-y plot

    # title = 'pitch = {}, q/m ={}'.format(pitchangle, m)
    # energy plot
    size = [12.8, 9.6]
    plt.figure(2, figsize=size)
    plt.plot(T_plot, ((Ke0-Ke)/Ke0))
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

    # 3d plot
    # python isnt suited for this, as it does not have a true 3d renderer

    '''
    def VtoT(T):  # output in ev
        Tev = T * 6.241509e18
        return Tev
    plt.figure(2,dpi = 300)
    plt.plot(T_plot, VtoT(Vx))
    plt.plot(T_plot, VtoT(Vy))
    plt.plot(T_plot, VtoT(Vz))
    plt.legend(['Ke_x', 'Ke_y', 'Ke_z'], loc = 'upper right')
    plt.title(filename)
    plt.xlabel(timelabel)
    plt.ylabel('Kinetic Energy (eV)')
    #plt.savefig('ParticlePlots/'\
        + filename+ 'CartesianVelocity.png' ,format = 'png')

    #V, Vparr,Vperp
    plt.figure(3,dpi = 300)
    #plt.plot(t,V0/V) #test for magnitude of V\
        after conversion to vparr and vperp
    plt.plot(T_plot, VtoT(Vparrallel))
    plt.plot(T_plot, VtoT(Vperpendicular))
    plt.legend([ 'T parallel', 'T perpendicular'], loc = 'upper right')
    plt.title(filename)
    plt.xlabel(timelabel)
    plt.ylabel('Energy (eV)')
    #plt.ylim(0,1.2)
    plt.savefig('ParticlePlots/'+filename + \
                'Parallel-Perpendicular Velocity.png' ,format = 'png')
'''


'''
    plt.figure(5,dpi = 300)
    (L,o) = trsfrm.cartesiantoLshell(xline,yline,zline)
    plt.plot(T_plot,L)
    #plt.plot(T_plot,o)
    plt.legend([ 'L '], loc = 'upper right')
    plt.title(filename)
    plt.xlabel(timelabel)
    plt.ylabel('L')
    #plt.savefig('ParticlePlots/'+filename+'L_shell.png' ,format = 'png')
    '''
'''
    ax = plt.axes(projection='3d')
    # generates earth size sphere,
    # looks weird though since perspective doesn't work
    u, v = np.mgrid[0:2*np.pi:200j, 0:np.pi:100j]
    x = np.cos(u)*np.sin(v)
    y = np.sin(u)*np.sin(v)
    z = np.cos(v)
    ax.plot_surface(x, y, z, color="r")

    ax.plot(xline, yline, zline, 'b')
    #saving 3d plots doesn't seem to work, So this will be left disabled'
'''

'''
def three_d_plot(xpoints,ypoints,zpoints,filename,title):
    # 3d animation, very resource intensive and limited length

    fig = plt.figure(figsize=(9, 9))

    axf = plt.subplot(3,3,(1,9),projection = '3d')

    # final graph
    axf.set_xlim(-2,2)
    axf.set_ylim(-2,2)
    axf.set_zlim(-1,1)
    axf.grid()
    axf.set_xlabel('x(m)')
    axf.set_ylabel('y(m)')
    axf.set_zlabel('z(m)')

    txt_title = axf.set_title('')

    linef, = axf.plot3D([],[],[],'-b',lw = 1)
    ptf, = axf.plot3D([],[],[],'.k', ms = 20)
    #tpoints = np.linspace(0,50,10000)
    #tpoints = t

    def drawframe(n):
        arraysize = len(xpoints)
        x = xpoints[n]
        y = ypoints[n]
        z = zpoints[n]
        #t = tpoints[n]
        txt_title.set_text(title)
        linef.set_data_3d(xpoints[0:n],ypoints[0:n],zpoints[0:n])
        ptf.set_data_3d(x,y,z)

        angle = (360/(arraysize*10000))*10*n+45
        axf.view_init(30, angle)
        return (linef,ptf)

    from matplotlib import animation

    anim = animation.FuncAnimation(fig, drawframe, frames = 600, interval=1,
                                   blit=True,  cache_frame_data = False)
    # run out of memory if frame number is too high > 1000?

    anim.save('ParticlePlots/'+filename+'_3ddipolemotion.mp4'
              , writer='imagemagick',fps = 24)

'''
