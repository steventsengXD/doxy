# Importing the necessary libraries 
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset


# Read in the fortran data 
filename='ising.nc'
f90data=Dataset(filename,mode='r')


# Retrieve and convert the simulation data 
N=f90data.N_gridsize  # The dimension of grid
t=np.array(f90data.variables['t'][:])  # t axis
netmag=np.array(f90data.variables['net_magnetization'][:])  # net magnetization
grid_final=np.array(f90data.variables['final_grid'][:])  # final grid 
grid_init=np.array(f90data.variables['initial_grid'][:])  # initial grid 


# Dictionary of initial configurations 
config={'A':'Alternating','a':'Alternating','F':'Ferromagnet',\
             'f':'ferromagnet','R':'Random','r':'Random'}


# Convert command line input to full word for plot title
grid0=config[f90data.init_config]


# Fontsizes for subplot titles and axes labels, respectively 
fs1=13  
fs2=11  


# Three subplots to visualize our simulation results 
# Some pretty colormaps we can swap in: Greys,cividis,Reds,bwr,Blues
fig,ax=plt.subplots(nrows=1,ncols=3,figsize=(18,5))

# Plot of initial grid
plt.subplot(1,3,1)
plt.imshow(grid_init, cmap='Blues',vmin=-1,vmax=1,extent=[1,N,N,1]) #,interpolation='gaussian') 
plt.title('{} Initial Grid, $t = 0$'.format(grid0),fontsize=fs1)
plt.xlabel(xlabel='Column $(i)$',fontsize=fs2)
plt.ylabel(ylabel='Row $(j)$',fontsize=fs2)
plt.colorbar()
plt.axis('tight')

# Plot of final grid
plt.subplot(1,3,2)
plt.imshow(grid_final,cmap='Blues',vmin=-1,vmax=1,extent=[1,N,N,1]) #,interpolation='gaussian') 
plt.title('Final Grid, $t = {}$'.format(len(netmag)-1),fontsize=fs2)
plt.xlabel(xlabel='Column $(i)$',fontsize=fs1)
plt.ylabel(ylabel='Row $(j)$',fontsize=fs1)
plt.colorbar()
plt.axis('tight')

# Plot of net magnetization
plt.subplot(1,3,3)
plt.xlabel(xlabel='timestep ($t$)',fontsize=fs1)
plt.ylabel(ylabel='Net Magnetization',fontsize=fs1)
plt.title('Ising Model Relaxation History',fontsize=fs2)
plt.plot(t,netmag)
plt.axis('tight')

plt.tight_layout()
plt.show()
