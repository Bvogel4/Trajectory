# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 11:17:51 2019

@author: Blake
"""

import numpy as np
import math as math
from matplotlib import pyplot as plt
from scipy.integrate import solve_ivp
from sklearn.preprocessing import normalize
from mpl_toolkits import mplot3d
import os
from scipy.fft import fft, fftfreq
from vtk_export import vtk_export
import tempfile
#creates new folder for animation if it doesn't already exist
if not os.path.exists('ParticlePlots'):
        
     os.mkdir('ParticlePlots')
#main functions 

    
''' algorithim loses accuracy at this scale 
Re = 6.37e6 # radius of earth in m

m = 9.11e-31 #mass
M = 8e22 #dipole moment?   M = 8x10^22 \A m^2 for earth?
e = 1.602e-19 #charge of particle

#intial conditions (can't be 0,0,0)
x0 = 2*Re
y0 = 0
z0 = 0
Vx0 = 0
Vy0 = 1
Vz0 = 0
'''


def B(x,y,z):
     #magnetic field of a dipole
     M = 100 #magnetic moment M = 8x10^22 \A m^2 for earth
     
     r = np.sqrt(np.power(x,2) + np.power(y,2) +np.power(z,2))
     
    #bx = 3M xz/r5  M is a constant
     Bx = 3*M*x*z/np.power(r,5)
     By = 3*M*y*z/np.power(r,5)
     Bz = M*(3*np.power(z,2) - np.power(r,2))/np.power(r,5)
    
     return Bx,By,Bz 
     #return 0,0,1 # test case with constant B vs known solution
     
     
     
#U = [x,y,z,Vx,Vy,Vz]
def dUdt(t,U): 
    [x,y,z,Vx,Vy,Vz] = U
    (Bx,By,Bz) = B(x,y,z)
    
    #dVx/dt = e/m*(VyBz - VzBy)
    dVxdt = e/m*(Vy*Bz - Vz*By)
    
    #dVy/dt = e/m *(VzBx - VxBz)
    dVydt = e/m *(Vz*Bx - Vx*Bz)
    
    #dVz/dt = e/m* (VxBy - VyBx)
    dVzdt = e/m* (Vx*By - Vy*Bx)
    U = [x,y,z,Vx,Vy,Vz]
    
    #=!how to find x,y,z?
    
    return U[3] , U[4], U[5], dVxdt, dVydt, dVzdt

#make sure V is not aligned with position only
def cartesiantoLshell(x,y,z):
   
    r = np.sqrt(np.power(x,2) + np.power(y,2) +np.power(z,2))
    r = r/Re #converts r to be in earth radii
    lambda0 = np.arctan2(z, np.sqrt(x**2 + y**2))
    denum = np.cos(lambda0)**2
    L = r/denum
    return L, lambda0

def Bnormal(Bx,By,Bz):
    BdotB = Bx*Bx+By*By+Bz*Bz
    #print(BdotB)
    f = np.sqrt(BdotB)      # for large r TypeError: loop of ufunc does not support argument 0 of type float which has no callable sqrt method
    
    Bx = Bx/f
    By = By/f
    Bz = Bz/f
    return Bx,By,Bz
    
def VparVperp(x,y,z,Vx,Vy,Vz): #retrns Vpar and Vperp from 6 inputs
    
    #compute the B vector at each point of the particle motion
    Bx,By,Bz = B(x,y,z) #B(x,y,z)
    
    #add test case to find if the number is too large
    
    #normalizng
    Bx,By,Bz = Bnormal(Bx, By, Bz)
 
    
    Vparrallel = Vx*Bx + Vy*By + Vz*Bz
    
    Vperpx = Vy*Bz - Vz*By
    Vperpy = Vz*Bx - Vx*Bz
    Vperpz = Vx*By - Vy*Bx
    V_perp = np.array([Vperpx,Vperpy,Vperpz])
    
    
    #V0 = np.sqrt(soln.y[3]**2+soln.y[4]**2+soln.y[5]**2)
    #Vperpendicular = np.linalg.norm(V_perp)
    Vperpendicular = np.sqrt(np.power(Vperpx,2) + np.power(Vperpy,2) + np.power(Vperpz,2))
    #print(Vperpendicular)
    #print(Vparrallel)
    
    #V = np.sqrt(Vparrallel**2+Vperpendicular**2)
    return Vparrallel, Vperpendicular ,V_perp # last value is a vector, first 2 are scalar magnitudes
     

def main(x0,y0,z0,Vx0,Vy0,Vz0):  
    
    
    # time length and accuracy ( scipydoes not care about dt too much)
    dt = d_t
    tfinal = t_runtime
    
    looplimit = int(tfinal/dt)
    looplimit = looplimit +1
    t = dt*np.linspace(0,looplimit -1,looplimit)
    t_span = (0.0, tfinal) # Provide solution over this time range
    y0 = [x0,y0,z0,Vx0,Vy0,Vz0]             # Initial conditions
     
    t_eval = dt*np.linspace(0, looplimit -1, looplimit)
    soln = solve_ivp( dUdt, t_span, y0, method='RK45', t_eval=t_eval)
    
    # kinetic energy over time should be constant, variances are due to error
    
    T = .5 * m *( np.power(soln.y[3],2) + np.power(soln.y[4],2) + np.power(soln.y[5],2) )
    T0 = T[0]

    #getting V parrallel and V perpendicular
    Vx = soln.y[3]
    Vy = soln.y[4]
    Vz = soln.y[5]
    
    xline = soln.y[0]
    yline = soln.y[1]
    zline = soln.y[2]

    (Vparrallel, Vperpendicular,V_perp) = VparVperp(xline,yline,zline,Vx, Vy, Vz) 
    
    
       
    #for vtk export
    points = np.column_stack([xline,yline,zline])
    tmpdir = tempfile.gettempdir()
    out_filename = os.path.join(tmpdir,'particle_motion.vtk')
    #out_filename = 'C:/Users/blake/.spyder-py3/.particle_motion.vtk'
    ftype = 'ASCII' 
    connectivity = {'LINES' : np.array([5,5,3])}

    vtk_export(out_filename, points,
                    dataset = 'POLYDATA',
                    connectivity = connectivity,
                    title='Title',
                    ftype=ftype,
                    debug = False)
    
    if plot == True:
        '''
        2D Plots:
        '''
        
        #energy plot
        plt.clf()
        plt.figure(1)
        plt.plot(t,T/T0)
        #plt.ylim(0,1.2)
        plt.legend(['Kinetic Energy'])
        plt.title('Energy vs Time')
        plt.xlabel('Time(s)')
        plt.ylabel('Ratio of T/T0')
        plt.savefig('ParticlePlots/Energy.png' ,format = 'png')
        
        #velocity plots
        #x,y,z
        plt.figure(2)
        plt.plot(t, soln.y[3])
        plt.plot(t, soln.y[4])
        plt.plot(t, soln.y[5])
        plt.legend(['Vx', 'Vy', 'Vz'])
        plt.title('Cartesian Velocity')
        plt.xlabel('Time(s)')
        plt.ylabel('Velocity (m/s)')
        plt.savefig('ParticlePlots/CartesianVelocity.png' ,format = 'png')
        
        #V, Vparr,Vperp
        plt.figure(3)
        #plt.plot(t,V0/V) #test for magnitude of V after conversion to vparr and vperp
        plt.plot(t, Vparrallel)
        plt.plot(t, Vperpendicular)
        plt.legend([ 'V parrallel', 'V perpendicular'])
        plt.title('Parrallel and Perpendicular Velocity')
        plt.xlabel('Time(s)')
        plt.ylabel('Velocity (m/s)')
        #plt.ylim(0,1.2)
        plt.savefig('ParticlePlots/Parallel-Perpendicular Velocity.png' ,format = 'png')
        
        #postion
        plt.figure(4)
        plt.plot(t, soln.y[0])
        plt.plot(t, soln.y[1])
        plt.plot(t, soln.y[2])
        plt.legend(['x','y','z'])
        plt.title('Position in Cartesian')
        plt.xlabel('Time(s)')
        plt.ylabel('Distance (m)')
        plt.savefig('ParticlePlots/Cartesian position.png' ,format = 'png')
        
        
        plt.figure(5)
        #L-shell
        (L,o) = cartesiantoLshell(soln.y[0],soln.y[1],soln.y[2])
        plt.plot(t,L)
        plt.legend([ 'L '])
        plt.title('L-Shell Coords')
        plt.xlabel('Time(s)')
        plt.ylabel('Radial Distance (m)')
        plt.savefig('ParticlePlots/L-shell.png' ,format = 'png')
        
        #xy plot
        plt.figure(6)
        plt.plot(soln.y[0],soln.y[1])
        plt.title('XY position')
        plt.xlabel('X (m)')
        plt.ylabel('Y (m)')
        plt.savefig('ParticlePlots/XvsY.png' ,format = 'png')
        plt.show()
        
    '''
    #3d plot
    ax = plt.axes(projection='3d')
    
    # Data for a three-dimensional line 
   
    #ax.plot3D(xline, yline, zline)
    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    x = 10*np.cos(u)*np.sin(v)
    y = 10*np.sin(u)*np.sin(v)
    z = 10*np.cos(v)
    ax.plot_surface(x, y, z, color="r", alpha = .3)
    plt.title('3D plot')
    '''
 
    
    
    '''
    #export to csv
    export = np.array( [[t_eval], [xline], [yline],[zline]])
    #print(export)
    export_reshaped = np.transpose(export.reshape(export.shape[0], -1))
    #print(export_reshaped)
    np.savetxt('file.csv',export_reshaped,delimiter='.')
    '''
    
    
    #animation maker
    '''
    start = 0
    stop = 500
    #add multithreading here
    for i in range(start,stop): #bottleneck in this loop
       
        xdata = xline[i]
        ydata = yline[i]
        zdata = zline[i]
        ax.scatter3D(xdata, ydata, zdata, c=zdata, marker ='.' );
              
        plt.savefig('3dplot/file_{0:03d}.png'.format(i),dpi = 600)
        plt.clf()
        ax = plt.axes(projection='3d')
        ax.plot_surface(x, y, z, color="r", alpha = .3)
        plt.title('3D plot')
    #gif maker 
    import imageio
    images = []
    for i in range(start,stop): #this one doesn't take long
        images.append(imageio.imread('3dplot/file_{0:03d}.png'.format(i)))
    imageio.mimsave('3dplot/movie.gif', images)
    '''
    '''
    # foruier anaylis of x
    
    #wx = fft(soln.y[0])
    
    # Number of samplepoints
    f_s = 1/dt
    t = t
    x = soln.y[0]
    
    fig, ax = plt.subplots()
    ax.plot(t, x)
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Signal amplitude');
    
    from scipy import fftpack
    
    X = fftpack.fft(x)
    freqs = fftpack.fftfreq(len(x)) * f_s
    
    fig, ax = plt.subplots()
    
    ax.loglog(freqs, np.abs(X),)
    ax.set_xlabel('Frequency in Hertz [Hz]')
    ax.set_ylabel('Frequency Domain (Spectrum) Magnitude')
    ax.set_xlim(0, 1)
    #ax.set_ylim(0, 10000)
    '''
  #trajectory- (x,y,z)=(%.1f,%.1f,%.1f); (vx,vy,vz)=(%.1f,%.1f,%.1f) GSM

#convert intital conditions to pitch angle, azimuth angle(phase?) magnetic moment and energy
#postion L shell and 2 angles

def car2ctd(x0,y0,z0,Vx0,Vy0,Vz0): 
    
    T = .5 * m *( np.power(Vx0,2) + np.power(Vy0,2) + np.power(Vz0,2) )
    
    #funtioncall to get Vpar and Vperp
    (Vparrallel, Vperpendicular, V_perp) = VparVperp(x0,y0,z0,Vx0,Vy0,Vz0)
    #pitch angle
    alpha = np.arctan2(Vperpendicular,Vparrallel)
    r = (x0,y0,z0)
    R_mag = np.sqrt(np.power(r[0],2) + np.power(r[1],2) + np.power(r[2],2)  )
    latitude = np.arcsin(z0/R_mag)
    longitude = np.arctan2(y0,x0)
    (L,latitude) = cartesiantoLshell(x0, y0, z0)
    (Bx , By, Bz) = B(x0,y0,z0)
    B0 = np.array([Bx , By, Bz])
    #phase shouldn't matter much for gyroid motion 
    rPar = (np.dot(r,B0)/(np.dot(B0,B0)))*B0
    rPerp = r - rPar
    #print(rPerp,'\n')
    rPerpM = np.sqrt(np.power(rPerp[0],2) + np.power(rPerp[1],2) + np.power(rPerp[2],2)  )
    
    #print(V_perp,'\n')
    temp1 = np.dot(rPerp,V_perp)
    temp2 = Vperpendicular * rPerpM
    
    phase = np.arccos(temp1/temp2) #!!!! this only returns the angle from 0 to pi --- should be fixed now
    #take cross prodcut of r_perp and Vperp to find sign of phase
    cross = np.cross(rPerp,V_perp)
    direction = np.sign(np.dot(B0,cross))
    if direction == -1:
        phase = -phase +2*np.pi #should now return phase from 0-2pi
         
    
    #print(phase)
 
    return T, math.degrees(alpha), math.degrees(phase), L, math.degrees(longitude), math.degrees(latitude)


def ctd2car(pitch, phase, Kinetic_energy, Lshell, latitude, longitude):
    r = Lshell* np.power(np.cos(latitude),2)
    print(r)
    phi = np.pi/2 - latitude
    
    x = r * np.sin(phi) * np.cos(longitude)
    y = r * np.sin(phi) * np.sin(longitude)
    z = r* np.cos(phi)
    (Bx, By, Bz) = B(x,y,z)
    (Bx, By, Bz) = Bnormal(Bx, By, Bz)
    
    
    V_mag = np.sqrt(2/m * Kinetic_energy) 
    V_par = np.cos(pitch)* V_mag
    V_perp = np.sin(pitch)* V_mag
    print(V_mag,V_par,V_perp)
    
    Vparx = Bx * V_par
    Vpary = By * V_par
    Vparz = Bz * V_par

    #phase shouldn't matter much for gyroid motion 
    
    B0 = np.array([Bx , By, Bz])
    
    r = (x,y,z)
    print(r)
    rPar = (np.dot(r,B0)/(np.dot(B0,B0)))*B0
    print(rPar)
    rPerp = r - rPar
    print(rPerp)
    rperp_mag = np.sqrt(np.power(rPerp[0],2) + np.power(rPerp[1],2) + np.power(rPerp[2],2)  ) 
    if rperp_mag == 0: #can't normalize a vector with length 0!
        rperpx = 0
        rperpy = 0
        rperpz = 0
    else:
        rperpx = rPerp[0]/rperp_mag
        rperpy = rPerp[1]/rperp_mag
        rperpz = rPerp[2]/rperp_mag
        
    r_Perp = np.array([rperpx,rperpy,rperpz])
    vhat = np.cos(phase)* r_Perp + np.sin(phase) *(np.cross(r_Perp,B0)) #normalized v perpendicular direction
    
    Vperpx = vhat[0]*V_perp
    Vperpy = vhat[1]*V_perp
    Vperpz = vhat[2]*V_perp
    ''' old code for hardcode phase
    Vperpx = rperpx*V_perp
    Vperpy = rperpy*V_perp
    Vperpz = rperpz*V_perp
    '''
    Vx = Vperpx + Vparx
    Vy = Vperpy + Vpary
    Vz = Vperpz + Vparz
    
    
    return x, y, z, Vx, Vy, Vz


def particle_demo(L_shell = 4 ,
                 latitude = 0, #all angles in degrees
                 longitude = 0 , 
                 pitch = 90 ,
                 phase = 0 ,
                 Kinetic_energy = 2 , # in eV/kg
                 mass = 1e-19 , #in Kg
                 charge = 5e-18 , # in coulumbs
                 t = 100, #time in seconds
                 dt = .001,
                 plots = True):
#def particle_demo(L_shell,latitude ,longitude ,pitch ,phase ,Kinetic_energy/m,charge,t,dt):
    #convert all angles to degrees
    latitude = math.radians(latitude)
    longitude = math.radians(longitude)
    pitch = math.radians(pitch)
    phase = math.radians(phase)
    
    global plot
    plot = plots
    
    global Re
    Re = 1
    global m 
    m = 1e-19 #mass
    global t_runtime
    t_runtime = t
    global d_t
    d_t = dt
    global e
   
    e = charge
    Kinetic_energy = Kinetic_energy * 1.60218e-19
    
    x,y,z,Vx,Vy,Vz = ctd2car(pitch, phase, Kinetic_energy, L_shell, latitude, longitude)
    
    main(x,y,z,Vx,Vy,Vz)




particle_demo()


'''
loops = 1000
t = np.linspace(0, 2*np.pi,loops)
x = np.sin(t)
y = np.cos(t)
phase = np.zeros(loops)
for i in range (0,loops):
    T , pitch, phase[i], L, lon, lat = cartesiantomagnocoords(20,0,0,x[i],y[i],0)
  
#print(phase)
plt.plot(t,phase)
        
    #print('Kinetic Energy = {0:.2f} \npitch agle = {1:.2f} \nPhase = {2:.2f}\nL = {3:.2f} \nlongitude = {4:.2f} \nLatitude = {5:.2f}'.format( T , pitch, phase, L, lon, lat  ))


'''



'''
inp = input("Enter A,B,C,D for different initial conditions \n")
if inp == 'A':
    #inputs (x0,y0,z0,Vx0,Vy0,Vz0)   
    print('Input accepted')
  
    main(20,0,0,0,1,0)
    T , pitch, phase, L, lon, lat = cartesiantomagnocoords(20,0,0,0,0,.1)
    print('Kinetic Energy = {0:.2f} \npitch agle = {1:.2f} \nPhase = {2:.2f}\nL = {3:.2f} \nlongitude = {4:.2f} \nLatitude = {5:.2f}'.format( T , pitch, phase, L, lon, lat  ))
    
elif inp == 'B':
    print('Input accepted')
    main(20,0,0,0,.1,.1)
    
elif inp == 'C':
    print('Input accepted')
    main(20,0,0,0,1,.1)
elif inp == 'D':
    print('Input accepted')
    main(20,0,10,0,1,0)
    
else :
    print('Invalid input')
'''