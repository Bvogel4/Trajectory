import numpy as np

import constants
from particle_sim import trajectory


parameters = constants.parameters
parameters.update({'L_shell':5})

accuracy = [1e2, 1e4, 1e2, 1e2, 1e4, 1e2]

method = ['boris', 'boris', 'rk45', 'boris', 'boris', 'rk45']
pitch_angle = [90, 90, 90, 10, 10, 10]
compute_efficiency = np.zeros(6)*np.nan

for a in range(len(pitch_angle)):
    parameters.update(
        {
            'method': method[a],
            'accuracy': accuracy[a],
            'pitch_anlge': pitch_angle[a]
        })

    compute_efficiency[a] = trajectory(parameters, 'method', plot=True)

print(compute_efficiency)