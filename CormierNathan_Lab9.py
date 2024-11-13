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
    coeff = t_step/(2*h) # coefficient used by lax method
    p_max = rho_sl[0]  # setting the max value of rho to 1st index of input initial condtions

    # Initializing density, flow, final output, x scale, and time scale arrays
    rho = np.zeros(shape=(N))
    flow = np.zeros(shape=(N))
    rho_final = np.empty((N,n_step+1))
    xplot = (np.arange(N)-1/2.)*h - L/2.         
    tplot = np.linspace(0,n_step*t_step,n_step+1)

    # Incorporating stoplight initial condtions for rho
    rho[int(N/4):int(N/2)] = p_max
    rho[0:int(N/4)] = rho_sl[1]
    rho[int(N/2):] = rho_sl[1]

    # Setting the initial state, rho(x,0)
    rho_final[:,0] = np.copy(rho)  

    # Periodic boundary conditions
    ip = np.arange(N) + 1  
    ip[N-1] = 0          # ip = i+1 with periodic b.c.
    im = np.arange(N) - 1  
    im[0] = N-1          # im = i-1 with periodic b.c.
    
    # Compute the flow = (Density)*(Velocity)
    flow[:] = rho[:] * (v_max*(1 - rho[:]/p_max))
    
    # Compute rho using the lax method
    rho[:] = .5*( rho[ip] + rho[im] ) - coeff*( flow[ip] - flow[im] )

    # Filling the output array with the calculated values of rho at all times
    rho = np.expand_dims(rho,axis=1)
    rho_final[:,1:n_step] = rho

    # return the solved values of rho over time
    return rho_final, xplot, tplot

# Part 2
# Show the output of your code (time evolution of the traffic density) both as a 2D contour
# plot, i.e. x on the horizontal axis and t on the vertical axis with ρ(x, t) as colour-coded
# contours (see Figure 7.20, page 186, in your text for an example), and as a “snapshot”
# plot in which in which ρ(x, t) is plotted versus x at selected time intervals, all on the
# same axes. Include these plots in your submission as
# LastnameFirstname_Lab9_Fig1.png and LastnameFirstname_Lab9_Fig2.png.

# calling function for desired input parameters
# and assigning the output values of rho, xplot, and tplot for graphing
r_xt, xp, tp = advection(divisions,length,vm,rho_y0)


print(np.flipud(np.rot90(r_xt)))
# 2D contour plotting code
# lvls = np.linspace(0., 1., num=11) 
# ct2d = plt.contour(xp, tp, np.flipud(np.rot90(r_xt)), lvls) 
# plt.clabel(ct2d, fmt='%1.2f') 
# plt.xlabel('x')
# plt.ylabel('time')
# plt.title('Density contours')
# plt.show()
