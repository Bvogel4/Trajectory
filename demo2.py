from particle_sim import trajectory

M_p = 1.6726219e-27  # kg
M_e = 9.10938356e-31  # kg
C_e = -1.60218e-19   # C


m = [M_e, M_e, M_p, M_p, M_p, M_p, M_p, M_p]
q = [C_e, C_e, -C_e, -C_e, -C_e, -C_e, -C_e, -C_e]
Ke = [1e3, 1e3, 1e8, 1e8, 1e8, 1e8, 1e8, 1e8]
# T = [1e-5, 1, 1e-2, .2, 10, 10, .1, .1]
T = [1e-5, 1, 1e-2, .2, 10, 10, .1, .1]
acc = [1e2, 1e2, 1e2, 1e3, 1e2, 1e2, 1e4, 1e2]
pitch = [90, 89, 90, 89, 90, 90, 5, 5]
Lshell = [2.1, 1, 1.64, 1, 3, 3, 1, 1]
sample = [20, .01, 10, 5, 1, 1, 1, 1]
inte = ['boris', 'boris', 'boris', 'boris', 'boris', 'rk45',
        'boris', 'rk45']

# takes my computer 30 mins to run all of these
# wow that last one takes ages
# why can't I parallize this? they share values so it comes out as nonsense

# Parallel(n_jobs=8, prefer="threads")(
#     delayed(trajectory)(pitch[i], m[i], Ke[i], q[i], T[i],
#                         acc[i], Lshell[i], sample[i], inte[i]) \
# for i in [0, 2, 3, 4, 5, 6])
# the short gyro and electron bounce period cause problems
# try reoganizing and doing them in serial while the others are parallel

for i in range(0, 1):
    trajectory(pitch[i], m[i], Ke[i], q[i], T[i],
               acc[i], Lshell[i], sample[i], inte[i])

    # len(m)):  # remeber to set this back to len(m)
