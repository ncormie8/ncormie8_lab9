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
    1500 numerical integration steps, along with
    appropriately sized x and t arrays for plotting
    the change in density with respect to position
    and time.'''

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
    rho_final = np.empty((N,n_step))
    xplot = np.arange(-N,N,h)        
    tplot = np.linspace(0,n_step*t_step,n_step)

    # Incorporating stoplight initial condtions for rho
    rho[int(N/4):int(N/2)] = p_max

    # Setting the initial state, rho(x,0)
    rho_final[:,0] = np.copy(rho)  

    # Periodic boundary conditions
    ip = np.arange(N) + 1  
    ip[N-1] = 0          # ip = i+1 with periodic b.c.
    im = np.arange(N) - 1  
    im[0] = N-1          # im = i-1 with periodic b.c.
    
    for istep in range(n_step):
     # Compute the flow = (Density)*(Velocity)
        flow[:] = rho[:] * (v_max*(1 - rho[:]/p_max))
    
    # Compute rho using the lax method
        rho[:] = .5*( rho[ip] + rho[im] ) - coeff*( flow[ip] - flow[im] )

    # Fill the output array with the calculated value of rho at time istep+1
        rho_final[:,istep] = np.copy(rho)
    
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
r_xt = np.flipud(np.rot90(r_xt))

# 2D contour plotting code
lvls = np.linspace(0., rho_m, num=11) 
ct2d = plt.contour(xp, tp, r_xt, lvls) 
plt.clabel(ct2d, fmt='%1.2f')
plt.colorbar()
plt.xlabel('x position')
plt.ylabel('time (s)')
plt.title('Density Contours')
plt.show()

# Snapshot plotting code

# Setting number of evenly place intervals / snapshots of rho to plot 
t_intervals = 5

# Initializing array of rho arrays to be plotted across the number of snaphots
rho_plot = np.zeros((t_intervals,600)) # ini

# Determining the step size between snapshots for plotting
stepsize = np.size(tp)/t_intervals

# Loop which divides tp into "t_intervals" evenly spaced time intervals of "stepsize" timesteps each
for p in range(t_intervals):  # chose p as looping variable since tp is being parsed
    rho_plot[p,:] = r_xt[(p*int(stepsize)),:]  # set rho to be plotted equal to rho at the evenly space time step index
    
    # Plot rho at time interval with conditional label against xp
    plt.plot(xp,rho_plot[p,:],label=('Interval '+ str(p+1)+'(t='+str(0.08*p*int(stepsize))+' s)')) 

# Conditionally generate title
plt.title('Time evolution of density over '+str(t_intervals)+' evenly spaced intervals')
# Place legend in best place for least overlap of graph data, set font small
plt.legend(loc=0,fontsize='small')
plt.xlabel('x position')
plt.ylabel('Density')
plt.show()

#link to my repository
#https://github.com/ncormie8/ncormie8_lab9