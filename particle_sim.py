import numpy as np
from datetime import datetime
from datetime import timedelta

import transformations as trsfrm
import integrator
import constants



def particle_sim(parameters = constants.parameters):
    
    
    if parameters['show_timing']:
        startTime = datetime.now()

    #print(L_shell,pitch_angle,mass,charge,t,Kinetic_energy,method,accuracy,
    #sampling,losscone,latitude,longitude,phase)
    # method,accuracy,sampling,losscone,show_timing = parameters['method'],\
    #     parameters['accuracy'],parameters['sampling'],\
    #     parameters['loss_cone'],parameters['show_timing']
    
    L_shell,mass,charge,t,species = parameters['L_shell'],parameters['mass'],\
        parameters['charge'],parameters['time'],parameters['species']
    
    # sampling modification
    sampling,accuracy = parameters['sampling'],parameters['accuracy']
    sampling /=  (np.pi*2) * L_shell**3
    accuracy/=  L_shell**3

    # internally all angles are in radians
    # latitude = np.radians(latitude)
    # longitude = np.radians(longitude)
    # phase = np.radians(phase)
    # pitch = np.radians(pitch_angle)

    # covert energy to joules
    #Kinetic_energy = Kinetic_energy * constants.C_e

    dT = 1/parameters['accuracy']

    # convert initial conditons to cartesian
    # x0, y0, z0, vx0, vy0, vz0 = trsfrm.ctd2car(pitch,
    #                                            phase, Kinetic_energy,
    #                                            L_shell, latitude, longitude,
    #                                            mass, constants.Re)
    x0, y0, z0, vx0, vy0, vz0 = trsfrm.ctd2car(parameters)
    S0 = np.array([x0, y0, z0, vx0, vy0, vz0])

    # Relativistic modifications
    v = np.linalg.norm([vx0, vy0, vz0])
    beta = v/constants.c
    gamma = (1-beta**2)**-.5
    mass = gamma*mass

    # convert into new dimensionless units
    # need constants first
    # Re is at the top
    tc = abs((mass/(charge*constants.Bc)))
    # in nd form can't have negative tc
    # qsign passes the sign of the charge
    # simply invert B if charge is negative.
    qsign = np.sign(charge)

    # convert to nd positions and velocities
    S0 = S0/constants.Re
    S0[3:6] = S0[3:6] * tc

    # compute dimensionless time
    tp = t/tc
    
    # check loss cone angle and exit if angle is in loss cone
    Ba = np.linalg.norm(trsfrm.B(0, 0, constants.za))
    if parameters['loss_cone'] is True:
        B = np.linalg.norm(trsfrm.B(S0[0], S0[1], S0[2]))
        term = np.sqrt(B/Ba)
        if abs(term) > 1 or L_shell < 1:
            print('particle will hit the atmosphere, skipping')
            return

        alphac = np.arcsin(term)
        if parameters['pitch'] < alphac:
            print('particle will hit the atmosphere, skipping')
            return

    # Integrate
    if parameters['method'] == 'rk45':
        x, y, z, vx, vy, vz, T = integrator.rk45_nd(dT, tp, S0, qsign)
    elif parameters['method'] == 'euler_cromer':
        x, y, z, vx, vy, vz, T = integrator.euler_cromer(dT, tp, S0, qsign)
    elif parameters['method'] == 'boris':
       # print(parameters)
        x, y, z, vx, vy, vz, T = integrator.boris(dT, parameters['sampling'],
                                                  S0, tp, qsign)
    else:
        print('invalid method')
        return

    if parameters['show_timing']:
        print('integration took {}'.format(datetime.now() - startTime))

    # convert time into seconds
    t = T*tc

    #convert velocites into Re/s
    vx = vx / tc
    vy = vy / tc
    vz = vz / tc

    return t, x, y, z, vx, vy, vz

import output

def trajectory(parameters,type = 'test',plot = False,save = False):
    
    startTime = datetime.now()    
    if type in ['test','method','trajectory']:
        #calculate estimated drift period       
        parameters.update({'time':1.1*trsfrm.t_d(parameters)})
    elif type == 'bounce':
        #calculate estimated bounce period 
        parameters.update({'time':20*trsfrm.t_b(parameters)})
    
    T, x, y, z, vx, vy, vz = particle_sim(parameters)
    
    if save:
        output.save(parameters,T, x, y, z, vx, vy, vz)
    if output:
        output.plot(parameters,T, x, y, z, vx, vy, vz)
        
    if type == 'test':
        return y[-1], T[-1]
    
    elif type == 'method':
        compute_time = timedelta.total_seconds(datetime.now() - startTime)
        Ke0 = np.sqrt( vx[0]**2 + vy[0]**2 + vz[0]**2)
        Kef = np.sqrt( vx[-1]**2 + vy[-1]**2 + vz[-1]**2)
        err_E = abs((Kef-Ke0)/Ke0)
        compute_efficiency = err_E/compute_time
        print(compute_efficiency)

        return compute_efficiency
