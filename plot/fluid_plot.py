#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import pandas as pd

inputfile = "../IO/input"
outputfile = "../IO/output"

inData = pd.read_csv(inputfile,delimiter= '=', names = ["variable", "value"])

xlength = inData[inData.variable.str.strip() == "xlength"].values[0,1]
ylength = inData[inData.variable.str.strip() == "ylength"].values[0,1]
imax = int(inData[inData.variable.str.strip() == "imax"].values[0,1])
jmax = int(inData[inData.variable.str.strip() == "jmax"].values[0,1])

#-------------------------
# Creates the mesh
dx = xlength/imax
dy = ylength/jmax
x = np.arange(0,xlength,dx)
# to show everything with pcolor plot -> add a last element
x = np.append(x,x[-1]+dx)
# Also simulation adds a cell before and a cell after:
x = np.insert(x,0,x[0]-dx)
x = np.append(x,x[-1]+dx)

# Same for y
y = np.arange(0,ylength,dy)
y = np.append(y,y[-1]+dy)
y = np.insert(y,0,y[0]-dy)
y = np.append(y,y[-1]+dy)

# Creates the meshgrid for the plot
X, Y = np.meshgrid(x,y)         # the y indexing is on dim 0 and x is on dim 1

#-------------------------
# Reshape output data (also need to delete the last element which is a new line and shows up as a nan)
infile = open(outputfile, "r")
lines = infile.read().strip().split("\n") # the .strip() eliminates blank lines at beginning or end of file, but watch out if a line is supposed to begin or end with whitespace like a tab or space

# Prepare figure
fig, ax = plt.subplots(1, 1)
ims = []

for line in lines:
    data = np.fromstring(line, dtype=float, sep=' ')
    data = np.reshape(data,(jmax+2,imax+2))
    # append a column and a row
    dataResh = np.zeros((jmax+3,imax+3))
    dataResh[:-1,:-1] = data

    # Plot time
    ims.append((plt.pcolor(X, Y, dataResh, cmap='RdBu'),))
    #c = ax.pcolor(X, Y, dataResh, cmap='RdBu')#, vmin=-5, vmax=10)


im_ani = animation.ArtistAnimation(fig, ims, interval=200, repeat_delay=3000,blit=True)
ax.set_title('pcolor')
plt.colorbar()
plt.show()

infile.close()
