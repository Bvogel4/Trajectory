import numpy as np
import matplotlib.pyplot as plt


from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

def B(x,y,z):
     #magnetic field of a dipole
     M = 100 #magnetic moment M = 8x10^22 \A m^2 for earth
     
     r = np.sqrt(np.power(x,2) + np.power(y,2) +np.power(z,2))
     
    #bx = 3M xz/r5  M is a constant
     Bx = 3*M*x*z/np.power(r,5)
     
     By = 3*M*y*z/np.power(r,5)
     Bz = M*(3*np.power(z,2) - np.power(r,2))/np.power(r,5)
    
     return Bx,By,Bz 
     #return np.array([0,0,1]) # test case with constant B vs known solution

dt = 1e-2;
mass = .1;
charge = 1.;
runtime = 100
duration = int(runtime/dt)



p = np.array([2., 0, 1.22e-16]);
v = np.array([5.66, 0, 0.]);

#B = np.array([0., 0., 1.]);
b = np.array([0.,0.,1.])
E = np.array([0., 0., 0.]);

S = np.zeros((duration,3)) 
V = np.zeros((duration,3)) 

for time in range(duration):
    
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


plt.plot(time,p[0])
#plt.plot(S[:,0],S[:,1],'k',linewidth=2.0); 
plt.figure(6)
plt.axes().set_aspect('equal', 'datalim')
plt.plot(S[:,0],S[:,1])
#plt.title('Initial pitch angle = {0:.2f} degrees'.format(pitchangle) )
plt.xlabel('X (m)')
plt.ylabel('Y (m)')

#plt.savefig('ParticlePlots/'+filename+'XvsY.png' ,format = 'png')        

plt.show()