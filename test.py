import numpy as np

from particle_sim import trajectory
import constants

#compares to a previous run to make sure nothing changes
#test is a proton and antiproton drift comparing the last x and t values


parameters = constants.parameters
parameters['L_shell'] = 4
test_values = np.zeros((2, 2))

for a in range(2):
    parameters['charge'] = constants.C_e * (-1)**(a)
    test_values[a] = trajectory(parameters)

test_values_i = np.loadtxt('test.txt')

#uncomment to rewrite intial values in test.txt
#np.savetxt('test.txt', test_values)

#assert test
np.testing.assert_almost_equal(test_values_i, test_values, decimal=13)


