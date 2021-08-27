import numpy as np
import matplotlib.pyplot as plt


from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
from numba import njit
def B(x,y,z):
     #magnetic field of a dipole
     M = 8e22 #magnetic moment M = 8x10^22 \A m^2 for earth
     
     r = np.sqrt(np.power(x,2) + np.power(y,2) +np.power(z,2))
     
    #bx = 3M xz/r5  M is a constant
     Bx = 3*M*x*z/np.power(r,5)
     
     By = 3*M*y*z/np.power(r,5)
     Bz = M*(3*np.power(z,2) - np.power(r,2))/np.power(r,5)
    
     return Bx,By,Bz 
     #return np.array([0,0,1]) # test case with constant B vs known solution
@njit
def B_nd_dipole(S): #spherica nd coords 
    r_bar,th,phi = S[0:3]
    B_dipole = np.array([np.cos(th),np.sin(th),0]) / r_bar**3
    return B_dipole

dt = 1e-4;
dT = 1e-2;
mass = .1;
charge = 1.;
runtime = 1
duration = int(runtime/dt)

Re = 6.37e6


p = np.array([1, 0.,0. ]);
v = np.array([0., 1.,0.]);

from Transformations import cartesian_to_spherical
p[0],p[1],p[2],v[0],v[1],v[2] = cartesian_to_spherical(p[0],p[1],p[2],v[0],v[1],v[2]) # important that p,v are not ints in this step
print(p,v)


#B = np.array([0., 0., 1.]);
b = np.array([0.,0.,1.])
E = np.array([0., 0., 0.]);

@njit
def loop(p,v):
    S = np.zeros((duration,3)) 
    V = np.zeros((duration,3)) 
    for time in range(duration):
        '''
        (b[0],b[1],b[2]) = B(p[0],p[1],p[2])
        
        t = charge / mass * b * 0.5 * dt;
        s = 2. * t / (1. + t*t);
        v_minus = v + charge / (mass ) * E * 0.5 * dt;
        v_prime = v_minus + np.cross(v_minus,t);
        v_plus = v_minus + np.cross(v_prime,s);
        v = v_plus + charge / (mass ) * E * 0.5 * dt;
        p += v * dt;
        S[time,:] = p;
        V[time,:] = v; 
        '''
    #nd and no E field
        b = B_nd_dipole(p)
        
        U =  b *dT;
        s = 2. * U / (1. + U*U);
        v_minus = v 
        v_prime = v_minus + np.cross(v_minus,U);
        v_plus = v_minus + np.cross(v_prime,s);
        v = v_plus 
        p = p+ v * dT;
        S[time,:] = p;
        V[time,:] = v; 
        
        return S,V
        
    
S,V = loop(p,v)
from Transformations import spherical_to_cartesian

S[:,0],S[:,1],S[:,2],V[:,0],V[:,1],V[:,2] = spherical_to_cartesian(S[:,0],S[:,1],S[:,2],V[:,0],V[:,1],V[:,2])
#plt.plot(time,p[0])
#plt.plot(S[:,0],S[:,1],'k',linewidth=2.0); 
plt.figure(6)
plt.axes().set_aspect('equal', 'datalim')
plt.plot(S[:,0],S[:,1])
#plt.title('Initial pitch angle = {0:.2f} degrees'.format(pitchangle) )
plt.xlabel('X (m)')
plt.ylabel('Y (m)')

#plt.savefig('ParticlePlots/'+filename+'XvsY.png' ,format = 'png')        

plt.show()