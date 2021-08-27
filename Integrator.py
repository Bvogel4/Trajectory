# -*- coding: utf-8 -*-
"""
Created on Wed Aug 11 20:40:22 2021

@author: blake
"""
import numpy as np
from numba import njit
from scipy.integrate import solve_ivp

import Transformations as trsfrm
from Transformations import B as B_nd_dipole
'''
def B(x,y,z): #magnetic field of a dipole or any other arbitary magnetic field
     M = 100 #magnetic moment M = 8x10^22 \A m^2 for earth
     
     r = np.sqrt(np.power(x,2) + np.power(y,2) +np.power(z,2))
     
    #bx = 3M xz/r5  M is a constant
     Bx = 3*M*x*z/np.power(r,5)
     
     By = 3*M*y*z/np.power(r,5)
     Bz = M*(3*np.power(z,2) - np.power(r,2))/np.power(r,5)
    
     return Bx,By,Bz 
     #return 0,0,1 # test case with constant B vs known solution
     
def B_nd_dipoleS(S): #spherical nd coords 
    r_bar,th,phi = S[0:3]
    B_dipole = np.array([2*np.cos(th),np.sin(th),0]) / r_bar**3

    #B[0],B[1],B[2]
    return B_dipole

def B_nd_dipole(P): #cartesian nd coords 
    r2 = (np.dot(P,P))
    B_dipole = np.array([ 3*P[0]*P[2] , 3*P[1]*P[2] , ( 3*P[2]-r2 ) ]) / np.power(r2,5/2)

    return B_dipole
'''
#
#note dimensionless position, velocity and t
#needs to be in cartesian before passing to integrator
def dUdt_nd(T,U):
    P = U[0:3] #position
    V = U[3:6] # velocity
    x,y,z = P
    
    B = B_nd_dipole(P[0],P[1],P[2]) # magnetic function (currently in cartesian)
    Bx,By,Bz = B 
    Vx,Vy,Vz = V
    DvDT = np.cross(V,B) # note cross only works in cartesian Does this work with the lists of values?
    '''
    #let me try explit defination of cross product. don't think that's the problem
    #dVx/dt = e/m*(VyBz - VzBy)
    dVxdt = (Vy*Bz - Vz*By)
    
    #dVy/dt = e/m *(VzBx - VxBz)
    dVydt = (Vz*Bx - Vx*Bz)
    
    #dVz/dt = e/m* (VxBy - VyBx)
    dVzdt =  (Vx*By - Vy*Bx)
    U = [x,y,z,Vx,Vy,Vz]
    '''

    return U[3] , U[4], U[5] , DvDT[0], DvDT[1], DvDT[2] # dVxdt, dVydt, dVzdt
#legacy
'''
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
    
    return U[3] , U[4], U[5], dVxdt, dVydt, dVzdt
'''


#old
'''
def rk45(dt,
         looplimit,
         tfinal,
         P):
        (x0,y0,z0,Vx0,Vy0,Vz0) = P
        t = dt*np.linspace(0,looplimit -1,looplimit)
        t_span = (0.0, tfinal) # Provide solution over this time range
        S0 = [x0,y0,z0,Vx0,Vy0,Vz0]             # Initial conditions
        
        t_eval = dt*np.linspace(0, looplimit -1, looplimit)
        soln = solve_ivp( dUdt, t_span, S0, method='RK45', t_eval=t_eval)
        
        xline = soln.y[0]
        yline = soln.y[1]
        zline = soln.y[2]
        Vx = soln.y[3]
        Vy = soln.y[4]
        Vz = soln.y[5]
        
        return xline,yline,zline,Vx,Vy,Vz,t
#nod dimensional function to solve in spherical
'''
def rk45_nd(dT,
            tfinal,
            S0):
    #I have an extra, uneeded parameter^
    
    

    T = dT*np.linspace(0,tfinal,int(tfinal/dT))#need to find appropiate dT and T ranges
    #maybe I sshould make new functions to transfrom to and from dimensionless form?
    
    T_span = (T[0],T[-1]) # Provide solution over this time range

    
    soln = solve_ivp( dUdt_nd, T_span, S0, method='RK45', t_eval=T)
    
    
    xline = soln.y[0]
    yline = soln.y[1]
    zline = soln.y[2]
    Vx = soln.y[3]
    Vy = soln.y[4]
    Vz = soln.y[5]

    
    #should do conersions after function? no I should do it before and after calling 
    '''
    u0 = 1.25663706212e-6
    #B0 = M * u0 / (4*pi*Re^3) M = dipole moment
    M = 100
    B0 = M * u0 / (4*np.pi*Re**3)
    #tau = Mass / ( q * B0)
    tau = m / (e * B0)
    t = tau * T
    '''
    return xline,yline,zline,Vx,Vy,Vz,T

#integrator boris method
#boris method cares very much what dt is
#boris has the benifit of conserving energy, and enables the computation with a electric field
#but is inaccuarte at low pitch angles (following the magnetic field lines)
#n = int(accuracy/20) #allows for more precise computation withough storing all the values/ saves ram
#tied to accuracy parameter so that points per tc remains constant at 20 points per loop 
#increase denominator for more points per loop
def boris(dT, #compute time step
          accuracy, #storage accuracy
          P,
          duration
          ):#initial conditions
#placeholder?

    
    (x0,y0,z0,Vx0,Vy0,Vz0) = P
    
    p = np.array([x0, y0,z0 ]);
    v = np.array([Vx0, Vy0,Vz0]);
    

    @njit
    def loop(p,v):
        j = 0
        n = int(accuracy/20)
        S = np.zeros((duration,3)) 
        V = np.zeros((duration,3)) 
        for time in range(duration):
        #nd and no E field
            b = B_nd_dipole(p)
            
            U =  b *dT;
            s = 2. * U / (1. + U*U);
            v_minus = v 
            v_prime = v_minus + np.cross(v_minus,U);
            v_plus = v_minus + np.cross(v_prime,s);
            v = v_plus 
            p = p+ v * dT;
            
            if np.mod(time,n) == 0:# this only grabs every nth value for storage
                S[j,:] = p; 
                V[j,:] = v; 
                j = j+1
                
            S[time,:] = p;
            V[time,:] = v; 
            
            
            return S,V
            
        
    S,V = loop(p,v)
    from Transformations import spherical_to_cartesian
    
    S[:,0],S[:,1],S[:,2],V[:,0],V[:,1],V[:,2] = spherical_to_cartesian(S[:,0],S[:,1],S[:,2],V[:,0],V[:,1],V[:,2])
#        return S,V
    
    S,V = loop(p,v)
    xline = S[:, 0]
    yline = S[:, 1]
    zline = S[:, 2]
    Vx = V[:,0]
    Vy = V[:,1]
    Vz = V[:,2]
    return xline,yline,zline,Vx,Vy,Vz

'''
    m =-1
    e = 1
    mass = m
    charge = e;
    m =-1
    e = 1
    
    n = int(accuracy/20)
    v = P[3:6]
    p = P[0:3]
    @njit
    def loop(p,v):
        
        
        lon0 = np.arctan2(p[1],p[0])
        #only want to compute 1 shell, so create conditions for comparison when 1 shell has been completed
        ni = int(np.floor(duration/n))
        S = np.zeros((ni,3)) 
        V = np.zeros((ni,3))
        #fill with nans so that unallocated values can easily be pruned
        #matplotlib automatically ignores nan values
        S[:] = np.nan
        V[:] = np.nan
        print(S.shape)
        #initialize B and E fields
        b = np.array([0.,0.,1.]) # doesnt really matter as it gets replaced by dipole later
        E = np.array([0., 0., 0.]);
        j = 0
        initial = True # variable that tracks when the particle has left the initial region
        overlap = .02 # overlap in the shell in rads
        for i in range(duration): 
            
            (b[0],b[1],b[2]) = B(p[0],p[1],p[2])
            
            #boris method calculation
            t = charge / mass * b * 0.5 * dt;
            s = 2. * t / (1. + t*t);
            v_minus = v + charge / (mass ) * E * 0.5 * dt;
            v_prime = v_minus + np.cross(v_minus,t);
            v_plus = v_minus + np.cross(v_prime,s);
            v = v_plus + charge / (mass ) * E * 0.5 * dt;
            p = p + v * dt;
            
            if np.mod(i,n) == 0:# this only grabs every nth value for storage
                S[j,:] = p; 
                V[j,:] = v; 
                j = j+1
            '''
#checks to make sure particle has left initial region
'''
            lon = np.arctan2(p[1],p[0])
            if initial == True:
                if lon >= 1/4*np.pi:
                    #print('condition')
                    initial = False
            #stops the loop if the particle has completed a shell
            if initial == False:
                dlon = abs(lon - lon0 - overlap)
                if dlon <= 1e-4:
                    #print(lon0,lon)
                    #print('break')
                    break
      
#        return S,V
    
    S,V = loop(p,v)
    xline = S[:, 0]
    yline = S[:, 1]
    zline = S[:, 2]
    Vx = V[:,0]
    Vy = V[:,1]
    Vz = V[:,2]
    return xline,yline,zline,Vx,Vy,Vz
'''