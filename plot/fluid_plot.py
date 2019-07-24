#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
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

# Prepare figure
fig, axs = plt.subplots(1, 1)

#-------------------------
# Reshape output data (also need to delete the last element which is a new line and shows up as a nan)
infile = open(outputfile, "r")
lines = infile.read().strip().split("\n") # the .strip() eliminates blank lines at beginning or end of file, but watch out if a line is supposed to begin or end with whitespace like a tab or space

for line in lines:
    data = np.fromstring(line, dtype=float, sep=' ')
    data = np.reshape(data,(imax+2,jmax+2))
    print(data.size)
    print(data.shape)
    # append a column and a row
    dataResh = np.zeros((imax+3,jmax+3))
    dataResh[:-1,:-1] = data

    # Plot time
    ax = axs
    c = ax.pcolor(X, Y, dataResh, cmap='RdBu')#, vmin=-5, vmax=10)
    ax.set_title('pcolor')
    fig.colorbar(c, ax=ax)
    plt.show()

infile.close()
