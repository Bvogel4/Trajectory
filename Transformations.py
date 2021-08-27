# -*- coding: utf-8 -*-
"""
Created on Wed Aug 11 20:35:34 2021

@author: blake
"""

#lowest level file, do not add fucntion calls from others files
import numpy as np
from numba import njit
import math
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
     '''
def B(x,y,z): 
    P = np.array([x,y,z])
    r2 = x*x + y*y + z*z
    B_dipole = np.array([ 3*P[0]*P[2] , 3*P[1]*P[2] , ( 3*P[2]-r2 ) ]) / np.power(r2,5/2)

    return B_dipole[0],B_dipole[1],B_dipole[2]
    #return 0,0,100000

@njit
def B_nd_dipole(S): #spherica nd coords 
    r_bar,th,phi = S[0:3]
    B_dipole = np.array[np.cos(th),np.sin(th),0] / r_bar**3
    return B_dipole

def cartesian_to_spherical(x, y, z, v_x, v_y, v_z): 
    """
    Utility function (jitted) to convert cartesian to spherical.
    This function should eventually result in Coordinate Transformation Graph!
    """
    hxy = np.hypot(x, y)
    r = np.hypot(hxy, z)
    theta = np.arctan2(hxy, z)
    phi = np.arctan2(y, x)
    
    n1 = x ** 2 + y ** 2
    n2 = n1 + z ** 2
    v_r = (x * v_x + y * v_y + z * v_z) / np.sqrt(n2)
    v_th = (z * (x * v_x + y * v_y) - n1 * v_z) / (n2 * np.sqrt(n1))
    v_p = -1 * (v_x * y - x * v_y) / n1

    return r, theta, phi, v_r, v_th, v_p


@njit
def spherical_to_cartesian(r, th, p, v_r, v_th, v_p):
    """
    Utility function (jitted) to convert spherical to cartesian.
    This function should eventually result in Coordinate Transformation Graph!
    """
    x = r * np.cos(p) * np.sin(th)
    y = r * np.sin(p) * np.sin(th)
    z = r * np.cos(th)
    v_x = (
        np.sin(th) * np.cos(p) * v_r
        - r * np.sin(th) * np.sin(p) * v_p
        + r * np.cos(th) * np.cos(p) * v_th
    )
    v_y = (
        np.sin(th) * np.sin(p) * v_r
        + r * np.cos(th) * np.sin(p) * v_th
        + r * np.sin(th) * np.cos(p) * v_p
    )
    v_z = np.cos(th) * v_r - r * np.sin(th) * v_th
    '''
    a_x = 0
    a_y = 0
    a_z = (a_r * np.cos(th) - v_r * np.sin(th) 
           - ( v_r * np.sin(th) + r * v_th**2 * np.cos(th) + r*np.sin(th) * a_th))
    '''

    return x, y, z, v_x, v_y, v_z


#!!!
# neeed to pass Re and B into this file
#*
def cartesiantoLshell(x,y,z,Re): # converts cartesian coords to L-shell
    r = np.sqrt(np.power(x,2) + np.power(y,2) +np.power(z,2))
    r = r/Re #converts r to be in earth radii : currently not in use
    lambda0 = np.arctan2(z, np.sqrt(x**2 + y**2))
    denum = np.cos(lambda0)**2
    L = r/denum
    return L, lambda0

def Bnormal(Bx,By,Bz): # normalizes magnetic field to find direction, used in other coord systems
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


#trajectory- (x,y,z)=(%.1f,%.1f,%.1f); (vx,vy,vz)=(%.1f,%.1f,%.1f) GSM

#converts Cartesion to more convientent to a more convienent system
#postion L shell and 2 angles
#@njit()
def car2ctd(x0,y0,z0,Vx0,Vy0,Vz0,m,Re): #doesnt work with vector lists.  m = mass
    
    T = .5 * m *( np.power(Vx0,2) + np.power(Vy0,2) + np.power(Vz0,2) )
    
    
    #funtioncall to get Vpar and Vperp
    (Vparrallel, Vperpendicular, V_perp) = VparVperp(x0,y0,z0,Vx0,Vy0,Vz0)
    #pitch angle
    alpha = np.arctan2(Vperpendicular,Vparrallel)
    r = (x0,y0,z0)
    R_mag = np.sqrt(np.power(r[0],2) + np.power(r[1],2) + np.power(r[2],2)  )
    latitude = np.arcsin(z0/R_mag)
    longitude = np.arctan2(y0,x0)
    (L,latitude) = cartesiantoLshell(x0, y0, z0,Re)
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
    
    phase = np.arccos(temp1/temp2) # this only returns the angle from 0 to pi --- should be fixed now
    #take cross prodcut of r_perp and Vperp to find sign of phase
    cross = np.cross(rPerp,V_perp)
    direction = np.sign(np.dot(B0,cross))
    if direction == -1:
        phase = -phase +2*np.pi #should now return phase from 0-2pi
         
    
    #print(phase)
 
    return T,  math.degrees(alpha), math.degrees(phase), L, math.degrees(longitude), math.degrees(latitude)

#converts back to cartesian
def ctd2car(pitch, phase, Kinetic_energy, Lshell, latitude, longitude,m,re):
    r = Lshell* np.power(np.cos(latitude),2) * re
    #print(r)
    phi = np.pi/2 - latitude
    
    x = r * np.sin(phi) * np.cos(longitude)
    y = r * np.sin(phi) * np.sin(longitude)
    z = r* np.cos(phi)
    (Bx, By, Bz) = B(x,y,z)
    (Bx, By, Bz) = Bnormal(Bx, By, Bz)
    
    
    V_mag = np.sqrt(2/m * Kinetic_energy) 
    V_par = np.cos(pitch)* V_mag
    V_perp = np.sin(pitch)* V_mag
    #print(V_mag,V_par,V_perp)
    
    Vparx = Bx * V_par
    Vpary = By * V_par
    Vparz = Bz * V_par

    #phase shouldn't matter much for gyroid motion 
    
    B0 = np.array([Bx , By, Bz])
    
    r = (x,y,z)
   # print(r)
    rPar = (np.dot(r,B0)/(np.dot(B0,B0)))*B0
   # print(rPar)
    rPerp = r - rPar
    #print(rPerp)
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
     