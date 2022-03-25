import numpy as np

import transformations as trsfrm
from particle_sim import particle_sim


#compares to a previous run to make sure nothing changes
#test is a proton and antiproton drift comparing the last x and t values
def demo_test(mass, charge):

    #calculate estimated bouce and drift period
    R = L_shell  # only valid when starting at equator
    pitch = np.radians(pitch_angle)
    # covert energy to joules
    Ke = Kinetic_energy * 1.602176565e-19
    phase, latitude, longitude = 0, 0, 0
    x0, y0, z0, vx0, vy0, vz0 = trsfrm.ctd2car(pitch,
                                               phase, Ke,
                                               L_shell, latitude, longitude,
                                               mass, Re)
    v = np.linalg.norm([vx0, vy0, vz0])
    beta = v/c
    if mass == M_p:
        Cd = 8.481
    elif mass == M_e:
        Cd = 1.557e4
    else:
        print('drift calculation only for electrons and protons')
    td = trsfrm.t_d(R, beta, np.radians(pitch_angle), Cd)
    T, xline, yline, zline, Vx, Vy, Vz = particle_sim(
        L_shell, pitch_angle, mass, charge, td, Kinetic_energy, method,
        accuracy, sampling, losscone=False)

    return yline[-1], T[-1]


Re = 6.371e6
M_p = 1.6726219e-27  # kg
M_e = 9.10938356e-31  # kg
C_e = -1.60218e-19   # C
c = 3e8

mass = [M_p, M_p]
charge = [-C_e, C_e]
Kinetic_energy = 1e8    # eV
pitch_angle = 89  # degress
L_shell = 4
method = 'boris'
accuracy = 1e+3
sampling = 36
test_values = np.zeros((2, 2))

for a in range(len(mass)):
    test_values[a] = demo_test(mass[a], charge[a])

test_values_i = np.loadtxt('test.txt')

#assert test
np.testing.assert_almost_equal(test_values_i, test_values, decimal=13)

#uncomment to rewrite intial values in test.txt
#np.savetxt('test.txt', test_values)
