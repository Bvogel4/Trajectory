from particle_sim import trajectory
import constants

parameters = constants.parameters
parameters.update({
                    'Kinetic_energy': 1e7,
                    'pitch_angle': 89,
                    'L_shell': 5,
                    'accuracy': 1e4
                   })

mass    = [constants.M_p, constants.M_e]
charge  = [-constants.C_e, constants.C_e]
species = ['proton', 'electron']

for a in range(len(mass)):
    parameters.update(
        {
            'mass': mass[a],
            'charge': charge[a],
            'species': species[a]
        })

    trajectory(parameters, 'bounce', plot=True)
