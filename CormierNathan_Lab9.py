import numpy as np
import matplotlib.pyplot as plt

# Part 1
# Adapt the traffic.py code into a function that solves the advection equation for the time
# evolution of ρ(x, t), using the Lax method with ρSL(x) as the initial condition. Integrate for
# at least 1500 time steps using the time step � = h/vmax.

# Defining constant parameters
vm = 25        # given v_max = 25
rho_m = 1      # given rho_max = 1

# Logical determination of N and L to get desired n_step = 1500
# want n_step = 1500
# v_max is constant at 25
# t_step = h/v_max = h/25 = (L/N)/25
# n_step = (L/4)/(v_max*tstep)
# = (L/4)/(25*(L/N)/25)
# = (L/4)/((25L/N)/25)
# = (L/4)/(L/N)
# 1500 = (L/4)/(L/N)
# 1500*(L/N) = L/4  (L's cancel)
# 1500/N = 1/4
# 6000/N = 1
# therefore N = 6000, L = any

def advection(grid_pts,sys_size,v_max,rho_sl):
    '''Solves the advection equation for the 
    time evolution of rho(x,t) using the lax 
    method with rho_sl as the initial condition.
    Returns function values of rho(x,t) over the
    calculated number of steps for integration.'''

    # Defining parameters
    N = grid_pts       # number of grid points
    L = sys_size       # system size (m)
    h = L/N            # grid spacing for BCs
    t_step = h/v_max   # setting timestep to recommended value
    n_step = (L/4)/(v_max*t_step)  # last car stops moving after n_step
    coeff = t_step/2*h # coefficient used by lax method

    


    return