from particle_sim import trajectory
import constants

parameters = constants.parameters
parameters.update(
                    {
                        'Kinetic_energy': 3e7,
                        'pitch_angle': 20,
                        'L_shell': 2,
                        'accuracy':1e3
                     }
                )

kwargs = {
            'traj_type': 'bounce',
            'compute': True,
            'save': True,
            'plot': True,
            'animation': True
        }

# After running set compute=False to speed things up
trajectory(parameters, **kwargs)
