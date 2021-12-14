from particle_sim import particle_sim
from output import save, plot

M_p = 1.6726219e-27  # kg
M_e = 9.10938356e-31 # kg
C_e = -1.60218e-19   # C

m = M_p
q = -C_e
Ke = 1e8    # eV
T = .2      # seconds
pitch = 90. # degress
Lshell = 3.

acc = 1e+4
sample = 5

methods = ['euler', 'boris', 'rk45']

for method in methods:
    print("Computing trajectory using " + method)

    t, xline, yline, zline, Vx, Vy, Vz, err_V = particle_sim(
        pitchangle=pitch, mass=m, Kinetic_energy=Ke,
        charge=q, t=T, accuracy=acc, L_shell=Lshell,
        sampling=sample, method=method)

    save(t, xline, yline, zline, Vx, Vy, Vz, Lshell, pitch, q,
         m, Ke, method)

    plot(Lshell, pitch, q, m, Ke, method, err_V, t=t,
         xline=xline, yline=yline, zline=zline,
         Vx=Vx, Vy=Vy, Vz=Vz)
