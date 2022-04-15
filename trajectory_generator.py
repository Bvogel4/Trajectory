import numpy as np

from joblib import Parallel
from joblib import delayed

import constants
from particle_sim import trajectory


def trajectory_generator(par=True):
    #start with quick ones in parrallel first 
    L_shell = np.linspace(10, 2, 5)
    pitch_angle = np.linspace(90, 15,6)
    mass = [constants.M_p, constants. M_e]
    charge = [-constants.C_e, constants.C_e]
    species = ['proton', 'electron']
    Kinetic_energy =  np.logspace(8, 3, 6)
    parameters = constants.parameters
    parameters.update({'losscone':True})

    if par is False:
        for a in range(len(mass)):
            for b in range(len(Kinetic_energy)):
                for c in range(len(pitch_angle)):
                    for d in range(len(L_shell)):
                        parameters.update({'L_shell': L_shell[d], 'pitch_angle':
                            pitch_angle[c], 'Kinetic_energy': Kinetic_energy[b],
                            'mass': mass[a], 'charge': charge[a],
                            'species': species[a]},plot = True,save = True)

                        trajectory(parameters)

    # use joblib to speed up compute time
    # saves both plots and trajectories, will take up a lot of space
    if par is True:
        Parallel(n_jobs=-2, prefer='threads')( #can't update dic like above
                                 #so have to input new dic
            delayed(trajectory)({'L_shell': 2, 'pitch_angle': 90, 'mass': 
                    mass[a], 'charge': charge[a], 'Kinetic_energy': 
                    Kinetic_energy[b], 'species': species[a], 'latitude': 0, 
                    'longitude': 0, 'phase': 0,'method': 'boris', 'accuracy':
                    1e3, 'sampling': 36, 'loss_cone': False, 'show_timing': 
                    False},save = True,plot = True)
            for a in range(len(mass))
            for b in range(len(Kinetic_energy))
            for c in range(len(pitch_angle))
            for d in range(len(L_shell)))


trajectory_generator(par=True)
