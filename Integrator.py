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

@njit
def dUdt_nd(T,U):
    #P = U[0],U[1],U[2] #position
    V = U[3],U[4],U[5] # velocity
    
    Bx,By,Bz = B_nd_dipole(U[0],U[1],U[2]) # magnetic field function (currently in cartesian)  numba points here
    B = Bx,By,Bz 
    Vx,Vy,Vz = V
    DvDT = np.cross(V,B) # note cross only works in cartesian 

    return U[3] , U[4], U[5] , DvDT[0], DvDT[1], DvDT[2] # dVxdt, dVydt, dVzdt

def rk45_nd(dT,
            tfinal,
            S0):
    #I have an extra, uneeded parameter^
    print('S0',S0)
    n = int(tfinal/dT)

    T = dT*np.linspace(0, n,n)#need to find appropiate dT and T ranges
    #maybe I sshould make new functions to transfrom to and from dimensionless form?
    
    T_span = (T[0],T[-1]) # Provide solution over this time range

    
    soln = solve_ivp( dUdt_nd, T_span, S0, method='RK45', t_eval=T,max_step = dT*100)
    
    
    xline = soln.y[0]
    yline = soln.y[1]
    zline = soln.y[2]
    Vx = soln.y[3]
    Vy = soln.y[4]
    Vz = soln.y[5]

    

    return xline,yline,zline,Vx,Vy,Vz,T

@njit
def euler_cromer(dT,
            tfinal,
            S0):
    
    n = int(tfinal/dT)
    w = np.zeros((n,6))
    
    w[0] = S0
    #print(w)
    T = np.linspace(0,n,n) * dT
    for j in range(0,n-1):
        dw =  np.array(dUdt_nd(dT,w[j]))  * dT
        w[j+1] = w[j] + dw 
 
    return w[:,0],w[:,1],w[:,2],w[:,3],w[:,4],w[:,5],T


#integrator boris method
#boris method with explicit dt
#boris has the benifit of conserving energy, and enables the computation with a electricfield - currently not enabled
#but is inaccuarte at low pitch angles (following the magnetic field lines)

def boris(dT, #compute time step
          sampling, #storage accuracy
          P, # intial conditions
          duration): # length of time to compute

    
    (x0,y0,z0,Vx0,Vy0,Vz0) = P
    
    p = np.array([x0, y0,z0 ]);
    v = np.array([Vx0, Vy0,Vz0]);
    

    @njit
    def loop(p,v):
        #allows for more precise computation withough storing all the values/ saves ram
        n = int(1/sampling/dT) # how many values to store base on sampling and dt
        #print(1/sampling /dT)
        loops = int(duration/dT)
        #without storing every point array size can be reduced
        arraysize = int(np.floor(loops/n) )
        S = np.zeros((arraysize,3)) 
        V = np.zeros((arraysize,3)) 
        T = np.zeros(arraysize)
        #fill with nans so that unused array slots aren't plotted
        T[:] = np.nan
        S[:] = np.nan
        V[:] = np.nan
        
        #T = np.linspace(0,loops,loops) * dT
        #t = np.linspace(0,n,duration) * dT * n
        j = 0
        k = 0
        initial = True
        lon0 = np.arctan2(p[1],p[0])
        overlap = .02 
        for time in range(loops):
        #nd and no E field
           
            Bx,By,Bz = B_nd_dipole(p[0],p[1],p[2])
            B = np.array([Bx,By,Bz])
            t = B * 0.5 * dT;
            s = 2. * t / (1. + t*t);
            v_minus = v 
            v_prime = v_minus + np.cross(v_minus,t);
            v_plus = v_minus + np.cross(v_prime,s);
            v = v_plus 
            p += v * dT;
            k = k + dT
            if np.mod(time,n) == 0:# this only grabs every nth value for storage
                S[j,:] = p; 
                V[j,:] = v; 
                T[j] = k
                j = j +1
                lon = np.arctan2(p[1],p[0])
                #code that stops copmutation after a full drift period
                #checks to see if particle has left initial region
                if initial == True:
                    if lon >= 1/4*np.pi:
                        #print('condition')
                        initial = False
                #stops the loop if the particle has completed a shell
                if initial == False:
                    dlon = abs(lon - lon0 - overlap)
                    if dlon <= 1e-4:
                        #print(lon0,lon)
                        print('break')
                        break
        return S,V,T
            
        
    S,V,T = loop(p,v)

    xline = S[:, 0]
    yline = S[:, 1]
    zline = S[:, 2]
    Vx = V[:,0]
    Vy = V[:,1]
    Vz = V[:,2]
    #print(T)
    
    return xline,yline,zline,Vx,Vy,Vz , T
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
'''