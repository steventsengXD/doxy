#################    PX913 MINI PROJECT      #################

#This part was written by Charlotte Rogerson
#For the PX913 Mini Project

#This code imports the NetCDF file and plots
#the Ex variable and plots the particle position
#in both x and y.


#Uploading the NetCDF file
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np


#Accessing the dataset
file = nc.Dataset('electron.nc4', mode='r', format="NETCDF4")
#Sanity checks - uncomment if not needed
print(file)

#Defining the variables in the file for the pseudocolour plot
#Only need Ex, Ey - Electric fields
Ex = file.variables['Ex_field']
Ey = file.variables['Ey_field']

#Sanity checks - uncomment if not needed
print(Ey)
print(Ex)

#Variables for the scatter plot
#Position
pos = file.variables['Position']

#Sanity checks - uncomment if not needed
print(pos)
#pos[:,0] = x_data, pos[:,1] = y_data,

#This sets the grid size - possibly change later to be user defined
N=Ex.shape[0]
xdata = np.linspace(-1, 1, N)
ydata = np.linspace(-1, 1, N)


#Pseudocolour plot of the variable Ex
#If it is being too slow use instead pcolormesh
#Need a 2d grid

fig,ax=plt.subplots(nrows=1,ncols=2,figsize=(12,5))

plt.subplot(1,2,1)
plt.pcolor(xdata, ydata, Ex, cmap = 'Blues' )
plt.title('Electric Field (Ex)')
plt.xlabel('$x$')
plt.ylabel('$y$')
plt.colorbar()
plt.axis('tight')


#Scatter graph of particle positioon in x vs y
plt.subplot(1,2,2)
plt.scatter(pos[0,1:], pos[1,1:])
plt.scatter(pos[1,0], pos[1,1], label='Initial Position')
plt.title('Electron Position')
plt.xlabel('$x$')
plt.ylabel('$y$')
plt.legend()
plt.axis('tight')

plt.tight_layout()
plt.show()



#
#
# plt.figure(1)
#
# plt.pcolor(xdata, ydata, Ex, cmap = 'Blues' )
#
# plt.title('Electric Field (Ex)')
# plt.xlabel('X')
# plt.ylabel('Y')
#
# plt.show()
#
#
# #Scatter graph of particle positioon in x vs y
# plt.figure(2)
#
# plt.scatter(pos[:,0], pos[:,1])
# plt.title('Electron Position')
# plt.xlabel('Electron position in X')
# plt.ylabel('Electron position in Y')
#
# plt.show()
