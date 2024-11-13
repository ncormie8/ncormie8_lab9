import numpy as np
import matplotlib.pyplot as plt

# Part 1
# Adapt the traffic.py code into a function that solves the advection equation for the time
# evolution of ρ(x, t), using the Lax method with ρSL(x) as the initial condition. Integrate for
# at least 1500 time steps using the time step � = h/vmax.

# Defining constant parameters
vm = 25        # given v_max = 25
rho_m = 1      # given rho_max = 1
length = 1200  # length of x grid from -600 to + 600 (given)
divisions = 600  # number of equally divided steps (given)
rho_y0 = [rho_m,0] # definition of the initial condition vector for rho(x,0)

def advection(grid_pts,sys_size,v_max,rho_sl):
    '''Solves the advection equation for the 
    time evolution of rho(x,t) using the lax 
    method with rho_sl as the initial condition.
    Returns function values of rho(x,t) over
    1500 numerical integration steps.'''

    # Defining parameters
    N = grid_pts       # number of grid points
    L = sys_size       # system size (m)
    h = L/N            # grid spacing for BCs
    t_step = h/v_max   # setting timestep to recommended value
    n_step = 1500      # minimum number of steps for integration
    coeff = t_step/2*h # coefficient used by lax method
    p_max = rho_sl[0]  # setting the max value of rho to 1st index of input initial condtions

    # Initializing density and flow arrays
    rho = np.zeros(shape=(N))
    flow = np.zeros(shape=(N))

    # Incorporating stoplight initial condtions
    rho[int(N/4):int(N/2)] = p_max
    rho[0:N/4],rho[N/2:] = rho_sl[1]

    # Periodic boundary conditions
    ip = np.arange(N) + 1  
    ip[N-1] = 0          # ip = i+1 with periodic b.c.
    im = np.arange(N) - 1  
    im[0] = N-1          # im = i-1 with periodic b.c.
    
    # Incorporating initial conditions





    return