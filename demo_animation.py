# -*- coding: utf-8 -*-
"""
Created on Fri Apr 15 22:51:23 2022

@author: blake
"""

from particle_sim import trajectory
import constants
from output import animation

parameters = constants.parameters
parameters.update({'Kinetic_energy':3e7,'pitch_angle':20,'L_shell':2,
                   'accuracy':1e3})

#after running set compute to false to speed things up
trajectory(parameters,traj_type='bounce',compute = True,save = True,
           plot = True,animation = True)
