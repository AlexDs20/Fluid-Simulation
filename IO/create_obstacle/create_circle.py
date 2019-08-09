#!/usr/bin/env python3.7
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import pandas as pd

inputfile = "../input"
outputfile = "../obstacle"

inData = pd.read_csv(inputfile,delimiter= '=', names = ["variable", "value"])

xlength = inData[inData.variable.str.strip() == "xlength"].values[0,1]
ylength = inData[inData.variable.str.strip() == "ylength"].values[0,1]
imax = int(inData[inData.variable.str.strip() == "imax"].values[0,1])
jmax = int(inData[inData.variable.str.strip() == "jmax"].values[0,1])

#-------------------------
# Creates the coordinates
dx = xlength/imax
dy = ylength/jmax
x = np.arange(0,xlength,dx)
y = np.arange(0,ylength,dy)

# Creates the mesh for the obstacles
# this only includes the inner domain, not the additional
# boundaries taht are added in the simulations

X, Y = np.meshgrid(x,y)         # Creates a mesh of coordinates
X = X + dx/2                    # Give the values of the center of the actual cells
Y = Y + dy/2

# IMPORTANT: Check the values to put in the include/definitions.h file
Obs = np.zeros(X.shape)

# For this example, we put an circle of radius 0.1 in the middle y and a quarter of the way along x
y0 = (y[0]+y[-1])*0.9
x0 = (x[0]+x[-1])*0.5
r = 0.05

circle = np.sqrt((X - x0) ** 2 + (Y - y0) ** 2)
circle = circle < r

# Optimally, we should read the definitions.h file here instead of hard coding the values
# Set the default obstacle value
Obs[circle] = 10

# An obstacle must be at least of 2 cell thickness to work.
# Hence, delete everything that is only one cell thick
for yi, yv in enumerate(y):
    for xi, xv in enumerate(x):
        if Obs[yi,xi]==10 and Obs[yi,xi-1]==0 and Obs[yi,xi+1]==0 or Obs[yi,xi]==10 and Obs[yi-1,xi]==0 and Obs[yi+1,xi]==0 :
            Obs[yi,xi] = 0

# Go throught the whole grid and assigne the more precise boundary values
# In principle going through the obstacles is enough but well...
for yi, yv in enumerate(y):
    for xi, xv in enumerate(x):
        if (Obs[yi,xi] == 10):
            if   ( Obs[yi+1,xi  ]==0 and  Obs[yi,xi-1]==0 ):   # NW boundary
                Obs[yi,xi] = 55
            elif ( Obs[yi+1,xi  ]==0 and  Obs[yi,xi+1]==0 ):   # NE boundary
                Obs[yi,xi] = 56
            elif ( Obs[yi-1,xi  ]==0 and  Obs[yi,xi-1]==0 ):   # SW boundary
                Obs[yi,xi] = 57
            elif ( Obs[yi-1,xi  ]==0 and  Obs[yi,xi+1]==0 ):   # SE boundary
                Obs[yi,xi] = 58
            elif ( Obs[yi+1,xi  ]==0 ):   # N boundary
                Obs[yi,xi] = 51
            elif ( Obs[yi-1,xi  ]==0 ):   # S boundary
                Obs[yi,xi] = 52
            elif ( Obs[yi  ,xi-1]==0 ):   # W boundary
                Obs[yi,xi] = 53
            elif ( Obs[yi  ,xi+1]==0 ):   # E boundary
                Obs[yi,xi] = 54

# Write the output file
file = open(outputfile,'w')
for yi, yv in enumerate(y):
    for xi, xv in enumerate(x):
        file.write( str(int(Obs[yi,xi])) + '\n' )
file.close()


# For good measure, show the plot
fig, ax = plt.subplots(1, 1)
h = ax.imshow(Obs, aspect='equal', origin='lower', interpolation='None', extent= (X[0,0], X[0,-1], Y[0][0], Y[-1,0]))
cb = fig.colorbar(h,ax=ax)
ax.title.set_text('Preview of obstacles')
ax.set_xlim(left=X[0,0],right=X[0,-1])
ax.set_ylim(bottom=Y[0,0],top=Y[-1,0])

for y_index, yv in enumerate(y):
    for x_index, xv in enumerate(x):
        if ( Obs[y_index, x_index]!=0 and Obs[y_index, x_index]!=10 ):
            label = int(Obs[y_index, x_index])
            text_x = xv + dx/2
            text_y = yv + dy/2
            ax.text(text_x, text_y, label, color='black', ha='center', va='center')

plt.show()
