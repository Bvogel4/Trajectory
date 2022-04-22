from particle_sim import trajectory
import constants

parameters = constants.parameters
parameters.update({'Kinetic_energy': 3e7, 'pitch_angle': 15, 'L_shell': 3,
                   'accuracy': 1e3})

kwargs = {
    'traj_type':'bounce',
    'compute':True,
    'save': True,
    'plot':False,
    'animation':True
    }

#after running set compute to false to speed things up
trajectory(parameters, **kwargs)
