# -*- coding: utf-8 -*-
"""
Created on Wed Aug 11 20:43:30 2021

@author: blake
"""
import numpy as np
from matplotlib import pyplot as plt
from vtk_export import vtk_export
import tempfile
import Transformations as trsfrm

def plotter(xline,yline,zline,Vx,Vy,Vz,filename,
            T, #nd time
            integrator,pitchangle,m,Re):
    Ke = .5 * m *( np.power(Vx,2) + np.power(Vy,2) + np.power(Vz,2) )
    Ke0 = Ke[0]
    


    #getting V parrallel and V perpendicular
    (Vparrallel, Vperpendicular,V_perp) = trsfrm.VparVperp(xline,yline,zline,Vx, Vy, Vz) 
    
    
    
    #for vtk export
    #xline[~np.isnan(xline) removes nan values, paraview can't work with nans
    #should probably just do this after collecting lists, so I don't have to worry about them later for some other function
    points = np.column_stack([xline[~np.isnan(xline)],yline[~np.isnan(yline)],zline[~np.isnan(zline)]])
    tmpdir = tempfile.gettempdir()
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
    ''' be done with integrator instead
    if integrator == 'boris':
      t = dt*np.linspace(0,looplimit -1,int(np.floor(looplimit/n)))
      
    Tc = 2 *np.pi * m / (e*Bi)
    timelabel = 'Time (Tc)'
    #more convient to plot in time scales involving gyro radius
    tc = t
    #tc = t/Tc
    '''
    timelabel = 'Time (Tc)'
    #energy plot
    plt.clf()
    plt.figure(1)
    if Ke0 == 0:
        print('divide by zero!')
        
    plt.plot(T,Ke/Ke0)
    #plt.ylim(0,1.2)
    plt.legend(['Kinetic Energy'])
    plt.title('Error in energy '+ filename+'method_' + integrator)
    plt.xlabel(timelabel)
    plt.ylabel('Ratio of T/T0')
    plt.savefig('ParticlePlots/' + filename + 'Energy.png' ,format = 'png')
    plt.ylim(.9, 1.1)
      
    #velocity plots
    #x,y,z
    def VtoT(T): #output in ev
        Tev = T * 6.241509e18
        return Tev
      
    plt.figure(2)
    plt.plot(T, VtoT(Vx))
    plt.plot(T, VtoT(Vy))
    plt.plot(T, VtoT(Vz))
    plt.legend(['Ke_x', 'Ke_y', 'Ke_z'])
    plt.title('Initial pitch angle = {0:.2f} degrees'.format(pitchangle) )
    plt.xlabel(timelabel)
    plt.ylabel('Kinetic Energy (eV)')
    plt.savefig('ParticlePlots/' + filename+ 'CartesianVelocity.png' ,format = 'png')
      
    #V, Vparr,Vperp
    plt.figure(3)
    #plt.plot(t,V0/V) #test for magnitude of V after conversion to vparr and vperp
    plt.plot(T, VtoT(Vparrallel))
    plt.plot(T, VtoT(Vperpendicular))
    plt.legend([ 'T parallel', 'T perpendicular'])
    plt.title('Initial pitch angle = {0:.2f} degrees'.format(pitchangle) )
    plt.xlabel(timelabel)
    plt.ylabel('Energy (eV)')
    #plt.ylim(0,1.2)
    plt.savefig('ParticlePlots/'+filename+'Parallel-Perpendicular Velocity.png' ,format = 'png')
      
    #postion
    plt.figure(4)
    plt.plot(T, xline)
    plt.plot(T, yline)
    plt.plot(T, zline)
    plt.legend(['x','y','z'])
    plt.title('Initial pitch angle = {0:.2f} degrees'.format(pitchangle) )
    plt.xlabel(timelabel)
    plt.ylabel('Distance (m)')
    plt.savefig('ParticlePlots/'+filename+'Cartesian position.png' ,format = 'png')
    plt.figure(5)
    
    #L-shell
    Re = 1
    (L,o) = trsfrm.cartesiantoLshell(xline,yline,zline)
    plt.plot(T,L)
    plt.legend([ 'L '])
    plt.title('Initial pitch angle = {0:.2f} degrees'.format(pitchangle) )
    plt.xlabel(timelabel)
    plt.ylabel('L')
    plt.savefig('ParticlePlots/'+filename+'L-shell.png' ,format = 'png')
    
    #xy plot
    plt.figure(6)
    plt.axes().set_aspect('equal', 'datalim')
    plt.plot(xline,yline)
    plt.title('Initial pitch angle = {0:.2f} degrees'.format(pitchangle) )
    plt.xlabel('X (m)')
    plt.ylabel('Y (m)')
    plt.savefig('ParticlePlots/'+filename+'XvsY.png' ,format = 'png')