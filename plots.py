# -*- coding: utf-8 -*-
"""
Created on Wed Aug 11 20:43:30 2021

@author: blake
"""
import numpy as np
from matplotlib import pyplot as plt
from vtk_export import vtk_export
import Transformations as trsfrm

def plotter(xline,yline,zline,Vx,Vy,Vz,filename,
            T,t, #nd and seconds time
            integrator,pitchangle,m,Re,units):
        #valide choices for units are 's' 'Tc' 'min' 'days'
    Ke = .5 * m *( np.power(Vx,2) + np.power(Vy,2) + np.power(Vz,2) )
    Ke0 = Ke[0]
    #plot in seconds
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
        print('invalid choice for units\nuse s,Tc,min, or days\ndefaulting to Tc')
        T_plot = T
        timelabel = 'Time (Tc)'
    #getting V parrallel and V perpendicular
    (Vparrallel, Vperpendicular,V_perp) = trsfrm.VparVperp(xline,yline,zline,Vx, Vy, Vz) 
    
    
    
    #for vtk export
    #xline[~np.isnan(xline) removes nan values, paraview can't work with nans
    #should probably just do this after collecting lists, so I don't have to worry about them later for some other function
    points = np.column_stack([xline[~np.isnan(xline)],yline[~np.isnan(yline)],zline[~np.isnan(zline)]])
    out_filename = (filename + 'particle_motion.vtk')
    #out_filename = 'C:/Users/blake/.spyder-py3/.particle_motion.vtk'
    ftype = 'ASCII'
    connectivity = {'LINES' : np.array([points.shape[0]])}

    vtk_export(out_filename, points, #should points be lines?
                    dataset = 'POLYDATA',
                    connectivity = connectivity,
                    title='Title',
                    ftype=ftype,
                    debug = False)
    #2d Plotting: Energies, position, L-shell,x-y plot

    #title = 'pitch = {}, q/m ={}'.format(pitchangle, m)
    #energy plot
    plt.figure(1,dpi = 300)
    plt.plot(T_plot,abs((Ke0-Ke)/Ke0))
    plt.title(filename)
    plt.xlabel(timelabel)
    plt.ylabel('Abslote error in energy')
    plt.savefig('ParticlePlots/' + filename + '_Energy.png' ,format = 'png')
    plt.show()
    plt.clf()
    #velocity plots
    #x,y,z
    def VtoT(T): #output in ev
        Tev = T * 6.241509e18
        return Tev
    '''
    plt.figure(2,dpi = 300)
    plt.plot(T_plot, VtoT(Vx))
    plt.plot(T_plot, VtoT(Vy))
    plt.plot(T_plot, VtoT(Vz))
    plt.legend(['Ke_x', 'Ke_y', 'Ke_z'], loc = 'upper right')
    plt.title(filename)
    plt.xlabel(timelabel)
    plt.ylabel('Kinetic Energy (eV)')
    #plt.savefig('ParticlePlots/' + filename+ 'CartesianVelocity.png' ,format = 'png')
      
    #V, Vparr,Vperp
    plt.figure(3,dpi = 300)
    #plt.plot(t,V0/V) #test for magnitude of V after conversion to vparr and vperp
    plt.plot(T_plot, VtoT(Vparrallel))
    plt.plot(T_plot, VtoT(Vperpendicular))
    plt.legend([ 'T parallel', 'T perpendicular'], loc = 'upper right')
    plt.title(filename)
    plt.xlabel(timelabel)
    plt.ylabel('Energy (eV)')
    #plt.ylim(0,1.2)
    plt.savefig('ParticlePlots/'+filename+'Parallel-Perpendicular Velocity.png' ,format = 'png')
      '''
    #position
    plt.figure(4,dpi = 300)
    plt.plot(T_plot, xline)
    plt.plot(T_plot, yline)
    plt.plot(T_plot, zline)
    plt.legend(['x','y','z'], loc = 'upper right')
    plt.title(filename)
    plt.xlabel(timelabel)
    plt.ylabel('Distance (Re)')
    plt.savefig('ParticlePlots/'+filename+'_Cartesian position.png' ,format = 'png')
    plt.show()
    plt.clf()
    #L-shell
    '''
    plt.figure(5,dpi = 300)
    (L,o) = trsfrm.cartesiantoLshell(xline,yline,zline)
    plt.plot(T_plot,L)
    #plt.plot(T_plot,o)
    plt.legend([ 'L '], loc = 'upper right')
    plt.title(filename)
    plt.xlabel(timelabel)
    plt.ylabel('L')
    #plt.savefig('ParticlePlots/'+filename+'L-shell.png' ,format = 'png')
    '''
    #xy plot
    plt.figure(6,dpi = 300)
    plt.axes().set_aspect('equal', 'datalim')
    plt.plot(xline,yline)
    plt.title(filename)
    plt.xlabel('X (Re)')
    plt.ylabel('Y (Re)')
    plt.savefig('ParticlePlots/'+filename+'_XvsY.png' ,format = 'png')
    plt.show()
    plt.clf()
    
    #3d plot 
    #python really isnt suited for this, as it does not have a true 3d renderer
'''
    ax = plt.axes(projection='3d')
    
    u, v = np.mgrid[0:2*np.pi:200j, 0:np.pi:100j]
    x = np.cos(u)*np.sin(v)
    y = np.sin(u)*np.sin(v)
    z = np.cos(v)
    ax.plot_surface(x, y, z, color="r")

    plt.axes(xline, yline, zline, 'b')
    
    zlim = np.max(xline)
    ax.set_zlim(-zlim, zlim)
    '''