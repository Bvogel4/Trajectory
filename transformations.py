# lowest level file, do not add function calls from other files
# contains useful functions such as B field, and coordinate conversions
import numpy as np
from numba import njit
import constants


@njit
def B(x, y, z):  # could add axis as direction of the dipole moment
    # non-dimensional dipole equation for magnetic field
    r2 = x**2 + y**2 + z**2
    d = np.power(r2, 5/2)

    Bx = 3 * x * z/d
    By = 3 * y * z/d
    Bz = (3 * z**2 - r2)/d

    return Bx, By, Bz
    # return B_dipole[0],B_dipole[1],B_dipole[2]
    # return 0,0,1


def cartesian_to_spherical(x, y, z, v_x, v_y, v_z):

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

    return x, y, z, v_x, v_y, v_z


# x,y,z should already be in units of re
def cartesiantoLshell(x, y, z):  # converts cartesian coords to L-shell
    r = np.sqrt(np.power(x, 2) + np.power(y, 2) + np.power(z, 2))
    # r = r/Re #converts r to be in earth radii
    lambda0 = np.arctan2(z, np.sqrt(x**2 + y**2))
    denum = np.cos(lambda0)**2
    L = r/denum
    return L, lambda0


def Bnormal(Bx, By, Bz):
    # normalizes magnetic field to find direction, used in other coord systems
    BdotB = Bx*Bx+By*By+Bz*Bz
    f = np.sqrt(BdotB)
    Bx = Bx/f
    By = By/f
    Bz = Bz/f
    return Bx, By, Bz


def VparVperp(x, y, z, Vx, Vy, Vz):  # retrns Vpar and Vperp from 6 inputs
    # compute the B vector at each point of the particle motion
    Bx, By, Bz = B(x, y, z)
    # normalizng
    Bx, By, Bz = Bnormal(Bx, By, Bz)
    Vparrallel = Vx*Bx + Vy*By + Vz*Bz
    Vperpx = Vy*Bz - Vz*By
    Vperpy = Vz*Bx - Vx*Bz
    Vperpz = Vx*By - Vy*Bx
    V_perp = np.array([Vperpx, Vperpy, Vperpz])
    # Vperpendicular = np.linalg.norm(V_perp)
    Vperpendicular = np.sqrt(np.power(Vperpx, 2) + np.power(Vperpy, 2)
                             + np.power(Vperpz, 2))
    # last value is a vector, first 2 are scalar magnitudes
    return Vparrallel, Vperpendicular, V_perp


# converts Cartesion to more convientent to a more convienent system
# postion L shell and 2 angles
# @njit()
def car2ctd(x0, y0, z0, Vx0, Vy0, Vz0, m, Re):
    # doesn't work with vector lists. m = mass
    T = .5 * m * (np.power(Vx0, 2) + np.power(Vy0, 2) + np.power(Vz0, 2))
    (Vparrallel, Vperpendicular, V_perp) = VparVperp(x0, y0, z0, Vx0, Vy0, Vz0)
    # pitch angle
    alpha = np.arctan2(Vperpendicular, Vparrallel)
    r = (x0, y0, z0)
    R_mag = np.sqrt(np.power(r[0], 2) + np.power(r[1], 2) + np.power(r[2], 2))
    latitude = np.arcsin(z0/R_mag)
    longitude = np.arctan2(y0, x0)
    (L, latitude) = cartesiantoLshell(x0, y0, z0, Re)
    (Bx, By, Bz) = B(x0, y0, z0)
    B0 = np.array([Bx, By, Bz])
    # phase shouldn't matter much for gyroid motion
    rPar = (np.dot(r, B0)/(np.dot(B0, B0))) * B0
    rPerp = r - rPar
    # print(rPerp,'\n')
    rPerpM = np.sqrt(np.power(rPerp[0], 2) + np.power(rPerp[1], 2)
                     + np.power(rPerp[2], 2))
    temp1 = np.dot(rPerp, V_perp)
    temp2 = Vperpendicular * rPerpM
    phase = np.arccos(temp1/temp2)
    # take cross prodcut of r_perp and Vperp to find sign of phase
    cross = np.cross(rPerp, V_perp)
    direction = np.sign(np.dot(B0, cross))
    if direction == -1:
        phase = -phase + 2 * np.pi

    return alpha, phase, T, L, latitude, longitude


# converts back to cartesian
def ctd2car(parameters):
    pitch, phase, Kinetic_energy, Lshell, latitude, longitude, m = \
        np.radians(parameters['pitch_angle']),np.radians(parameters['phase']),\
        parameters['Kinetic_energy']*abs(constants.C_e),parameters['L_shell'],\
        np.radians(parameters['latitude']),np.radians(parameters['longitude']),\
        parameters['mass']
        
            
    #print(Lshell,latitude,Re)
    r = Lshell * np.power(np.cos(latitude), 2) * constants.Re
    phi = np.pi/2 - latitude
    x = r * np.sin(phi) * np.cos(longitude)
    y = r * np.sin(phi) * np.sin(longitude)
    z = r * np.cos(phi)
    Bx, By, Bz = B(x, y, z)
    Bx, By, Bz = Bnormal(Bx, By, Bz)
    #c = 3e8
    sqrt = np.sqrt(2*constants.c**2*m + Kinetic_energy)

    #print(Kinetic_energy,constants.c,m)
    V_mag = np.sqrt(Kinetic_energy) * constants.c * sqrt / (
        constants.c**2*m + Kinetic_energy)

    V_par = np.cos(pitch) * V_mag
    V_perp = np.sin(pitch) * V_mag

    Vparx = Bx * V_par
    Vpary = By * V_par
    Vparz = Bz * V_par
    B0 = np.array([Bx, By, Bz])
    r = (x, y, z)
    rPar = (np.dot(r, B0)
            / (np.dot(B0, B0))) * B0
    rPerp = r - rPar
    rperp_mag = np.sqrt(np.power(rPerp[0], 2) + np.power(rPerp[1], 2)
                        + np.power(rPerp[2], 2))
    if rperp_mag == 0:  # can't normalize a vector with length 0!
        rperpx = 0
        rperpy = 0
        rperpz = 0
    else:
        rperpx = rPerp[0]/rperp_mag
        rperpy = rPerp[1]/rperp_mag
        rperpz = rPerp[2]/rperp_mag

    r_Perp = np.array([rperpx, rperpy, rperpz])
    # normalized v perpendicular direction
    vhat = np.cos(phase) * r_Perp + np.sin(phase) * (np.cross(r_Perp, B0))

    Vperpx = vhat[0]*V_perp
    Vperpy = vhat[1]*V_perp
    Vperpz = vhat[2]*V_perp

    Vx = Vperpx + Vparx
    Vy = Vperpy + Vpary
    Vz = Vperpz + Vparz

    return x, y, z, Vx, Vy, Vz

#equatorial bounce period
#low pitch angles will differ
def t_b(parameters):  # good to .5%
    #R0 is the distance from origin to the equatorial crossing
    #R = R0/Re
    #beta = v/c
    
    R = parameters['L_shell'] # only valid when starting at equator
    x0, y0, z0, vx0, vy0, vz0 = ctd2car(parameters)
    v = np.linalg.norm([vx0, vy0, vz0])
    beta = v/constants.c
    a_eq = np.radians(parameters['pitch_angle'])
    t_b = .117*R/beta * (1-.4635 * (np.sin(a_eq))**(3/4))
    return t_b

#equatorial drift period
#low pitch angles will differ
def t_d(parameters):
    
    R = parameters['L_shell'] # only valid when starting at equator
    x0, y0, z0, vx0, vy0, vz0 = ctd2car(parameters)
    v = np.linalg.norm([vx0, vy0, vz0])
    beta = v/constants.c
    
    
    if parameters['species'] == 'proton':
        Cd = 8.481
    elif parameters['species'] == 'electron':
        Cd = 1.557e4
    else:
        print('drift calculation only for electrons and protons')
        return
    
    gamma = (1-beta**2)**(-.5)
    #cd for e- = 1.557e4s
    # for proton = 8.481s
    a_eq = np.radians(parameters['pitch_angle'])
    t_d = Cd/R / (gamma*beta**2) * (1-.333 * np.sin(a_eq)**(.62))
    return t_d
