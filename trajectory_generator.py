import numpy as np
from datetime import datetime

from joblib import Parallel
from joblib import delayed

import constants
import transformations as trsfrm
import integrator
from output import save
from output import plot
from particle_sim import trajectory


    
    
    

def trajectory_generator(par=True):
    L_shell = np.linspace(2, 10, 2)
    # pitch = np.linspace(90, 10, 8)
    pitch_angle = [90]
    mass = [constants.M_p,constants. M_e]
    charge = [-constants.C_e, constants.C_e]
    species = ['proton','electron']
    # K = np.logspace(8, 2, 7)
    Kinetic_energy = [1e8]
    parameters = constants.parameters

    
    if par is False:
        for a in range(0,1): # len(mass)):  # skip electron for now
            for b in range(0, len(Kinetic_energy)):
                for c in range(0, len(pitch_angle)):
                    for d in range(0, len(L_shell)):
                        parameters.update({'L_shell':L_shell[d],'pitch_angle':
                           pitch_angle[c],'Kinetic_energy':Kinetic_energy[b],
                           'mass':mass[a],'charge':charge[a],'species':species[a]})
                        trajectory(parameters)
                        


    # use joblib to speed up compute time
    # carful when running this, I think if memory gets cached, it will break
    # energy error plots
    if par is True:
        Parallel(n_jobs=-2, prefer='threads')(
            delayed(trajectory)({'L_shell':2,'pitch_angle' : 90,'mass' : mass[a], 
                    'charge' : charge[a],'Kinetic_energy':Kinetic_energy[b],
                    'species': species[a],'latitude':0,'longitude':0,'phase':0,
                    'method' :'boris', 'accuracy' : 1e3,'sampling':36,
                    'loss_cone':False,'show_timing':False}) 
                                for a in range(0, len(mass)) 
                                for b in range(len(Kinetic_energy))
                                for c in range(len(pitch_angle)) 
                                for d in range(len(L_shell))) 
        
trajectory_generator(par = True)
