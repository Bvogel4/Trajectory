
import numpy as np
from numba import njit
from scipy.integrate import solve_ivp
import multiprocessing
import psutil
from transformations import B as B_nd_dipole

#cpu and memory management
n_threads = multiprocessing.cpu_count()
mem_t = psutil.virtual_memory().total


@njit
def dUdt_nd(T, U):
    V = U[3], U[4], U[5]  # velocity

    Bx, By, Bz = B_nd_dipole(U[0], U[1], U[2])
    # magnetic field function (currently in cartesian)
    B = Bx, By, Bz
    # Vx, Vy, Vz = V
    DvDT = np.cross(V, B)  # note cross only works in cartesian

    return U[3], U[4], U[5], DvDT[0], DvDT[1], DvDT[2]  # dVxdt, dVydt, dVzdt


@njit
def dUdt_nd2(T, U):
    # copy of the first function but negative for negative charges
    V = U[3], U[4], U[5]  # velocity

    Bx, By, Bz = B_nd_dipole(U[0], U[1], U[2])
    B = (-Bx, -By, -Bz)

    # Vx,Vy,Vz = V
    DvDT = np.cross(V, B)

    return U[3], U[4], U[5], DvDT[0], DvDT[1], DvDT[2]  # dVxdt, dVydt, dVzdt


def rk45_nd(dT, tfinal, S0, qsign):
    # create a t array for solution
    n = int(tfinal/dT)
    T = dT * np.linspace(0, n, n)  # need to find appropiate dT and T ranges

    #need to cap arraysize
    arraysize = len(T)
    #if 7*arraysize*8 >1e9
    T_span = (T[0], T[-1])  # Provide solution over this time range
    if qsign == 1:
        soln = solve_ivp(dUdt_nd, T_span, S0, method='RK45', t_eval=T,
                         atol=1e-13, rtol=1e-13)
    elif qsign == -1:
        soln = solve_ivp(dUdt_nd2, T_span, S0, method='RK45', t_eval=T,
                         atol=1e-13, rtol=1e-13)

    xline = soln.y[0]
    yline = soln.y[1]
    zline = soln.y[2]
    Vx = soln.y[3]
    Vy = soln.y[4]
    Vz = soln.y[5]

    return xline, yline, zline, Vx, Vy, Vz, T


@njit
def euler_cromer(dT, tfinal, S0, qsign):
    n = int(tfinal/dT)
    w = np.zeros((n, 6))
    w[0] = S0
    T = np.linspace(0, n, n) * dT

    if qsign == 1:
        f = dUdt_nd
    elif qsign == -1:
        f = dUdt_nd2

    for j in range(0, n-1):
        dw = np.array(f(dT, w[j])) * dT
        w[j+1] = w[j] + dw

    return w[:, 0], w[:, 1], w[:, 2], w[:, 3], w[:, 4], w[:, 5], T


'''
integrator boris method
boris method with explicit dt
boris has the benifit of conserving energy,
and enables the computation with a electric field - currently not enabled
but is inaccuarte at low pitch angles (following the magnetic field lines)
becasue of no adaptive step
'''


def boris(dT, sampling, P, duration, qsign):

    (x0, y0, z0, Vx0, Vy0, Vz0) = P

    p = np.array([x0, y0, z0])
    v = np.array([Vx0, Vy0, Vz0])

    @njit(nogil=True)
    def loop(p, v, dt):
        dT = dt
        # allows for more precise computation without storing all the values
        n = int((1/sampling/dT))
        # how many values to store base on sampling and dt
        loops = int(duration/dT)
        # without storing every point array size can be reduced
        arraysize = int(np.floor(loops/n))
        '''
        change this value to increase length of the array
        note the total memory used will be this value *7 since
        there are two 2darrays with added axis of length 3
        and another one for time
        #note that parallel computing requires additional ram for each job
        '''

        maxarraysize = mem_t/n_threads/1e9 / 7  # GB 7 arrays in total

        store_type = np.float64
        # changing this to np.float 32 is noticable in electron gyros

        # if changing number type change this too
        numbersize = 8  # itemsize does not work with numba?

        if (arraysize*7*numbersize) > (maxarraysize*1e9):
            # truncate time to keep accuracy of gyro
            arraysize = int(maxarraysize*1e9/numbersize)
            loops = int(arraysize*n)
            print('array has been truncated at size', 8*arraysize/1e9*7, 'GB')
        if loops > 6.6e9:
            loops = 6*10**9

        S = np.zeros((arraysize, 3), dtype=store_type)
        S[:] = np.nan
        V = np.copy(S)
        T = np.zeros(arraysize, dtype=store_type)
        T[:] = np.nan
        j = 0  # indexer
        k = 0  # time
        initial = True
        lon0 = np.arctan2(p[1], p[0])
        overlap = 1e-4
        temp_n = n

        for time in (range(loops)):
            # nd and no E field
            Bx, By, Bz = B_nd_dipole(p[0], p[1], p[2])
            B = np.array([Bx, By, Bz])
            if qsign == -1:
                B = -B
            t = B * 0.5 * dT
            s = 2. * t / (1. + t*t)
            v_minus = v
            v_prime = v_minus + np.cross(v_minus, t)
            v_plus = v_minus + np.cross(v_prime, s)
            v = v_plus
            p += v * dT
            k = k + dT
            if np.mod(time, temp_n) == 0:  # only grabs every nth # for storage
                # modifying dt and n to increase accuracy based on latitude
                R = (p[1]**2 + p[2]**2)**(.5)
                lat = np.arcsin(p[2]/R)
                f = abs(int(np.pi/2/(np.abs(lat) - np.pi/2)))
                #f = 1000*np.exp(-1*np.abs( (lat- np.pi/2)/5)) + 1
                temp_n = n * f
                dT = dt/f

                S[j, :] = p
                V[j, :] = v
                T[j] = k
                j = j + 1
                lon = np.arctan2(p[1], p[0])
                # code that stops copmutation after a full drift period
                # checks to see if particle has left initial region
                #x, y, z = p[0], p[1], p[2]
                #theta = abs((np.arctan2(np.sqrt(x**2+y**2), z)))
                #factor = np.sin(theta)**2
                #dT = dt * factor

                if initial is True:
                    if abs(lon0 - lon) >= np.pi/2:
                        initial = False
                # stops the loop if the particle has completed a shell
                if initial is False:
                    dlon = abs(lon - lon0 - overlap)
                    if dlon <= 1e-2:
                        # print(lon0,lon)
                        #print('1 drift period completed.
                        # Stopping simulation.')
                        break
        return S, V, T

    S, V, T = loop(p, v, dT)
    #
    # seperate coords
    xline = (S[:, 0])
    yline = (S[:, 1])
    zline = (S[:, 2])
    Vx = (V[:, 0])
    Vy = (V[:, 1])
    Vz = (V[:, 2])

    nans = np.isnan(xline)

    xline = xline[~nans]
    yline = yline[~nans]
    zline = zline[~nans]
    Vx = Vx[~nans]
    Vy = Vy[~nans]
    Vz = Vz[~nans]
    T = T[~nans]

    return xline, yline, zline, Vx, Vy, Vz, T
