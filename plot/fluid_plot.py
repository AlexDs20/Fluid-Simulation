#!/usr/bin/env python3.7
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
Xpc, Ypc = np.meshgrid(x,y)         # the y indexing is on dim 0 and x is on dim 1
Xq = Xpc + dx/2
Yq = Ypc + dy/2

#-------------------------
# Reshape output data (also need to delete the last element which is a new line and shows up as a nan)

Ufile = open(outputfile+'U', "r")
Ulines = Ufile.read().strip().split("\n") # the .strip() eliminates blank lines at beginning or end of file, but watch out if a line is supposed to begin or end with whitespace like a tab or space
Vfile = open(outputfile+'V', "r")
Vlines = Vfile.read().strip().split("\n")
Pfile = open(outputfile+'P', "r")
Plines = Pfile.read().strip().split("\n")

# Prepare figure
fig, ax = plt.subplots(2, 2)

for i,line in enumerate(Ulines):
    Udata = np.fromstring(line, dtype=float, sep=' ')
    Udata = np.reshape(Udata,(jmax+2,imax+2))
    # append a column and a row
    UdataResh = np.zeros((jmax+3,imax+3))
    UdataResh[:-1,:-1] = Udata

    Vdata = np.fromstring(Vlines[i], dtype=float, sep=' ')
    Vdata = np.reshape(Vdata,(jmax+2,imax+2))
    VdataResh = np.zeros((jmax+3,imax+3))
    VdataResh[:-1,:-1] = Vdata

    Pdata = np.fromstring(Plines[i], dtype=float, sep=' ')
    Pdata = np.reshape(Pdata,(jmax+2,imax+2))
    PdataResh = np.zeros((jmax+3,imax+3))
    PdataResh[:-1,:-1] = Pdata

    # Plot
#       if (i==0):
#           q = ax.quiver(Xq,Yq,UdataResh,VdataResh)
#       else:
#           q.set_UVC(UdataResh,VdataResh)
#       plt.draw()
    if (i==0):
        imU = ax[0,0].imshow(UdataResh, vmin=np.amin(Udata), vmax=np.amax(Udata), aspect='equal', origin='lower', extent= (Xpc[0,0], Xpc[0,-1]+dx, Ypc[0][0], Ypc[-1,0]+dy))
        cb = fig.colorbar(imU,ax=ax[0,0])
        ax[0,0].title.set_text('U')
        ax[0,0].set_xlim(left=Xpc[0,0],right=Xpc[0,-1])
        ax[0,0].set_ylim(bottom=Ypc[0,0],top=Ypc[-1,0])

        imV = ax[1,0].imshow(VdataResh, vmin=np.amin(Vdata), vmax=np.amax(Vdata), aspect='equal', origin='lower', extent= (Xpc[0,0], Xpc[0,-1]+dx, Ypc[0][0], Ypc[-1,0]+dy))
        cb = fig.colorbar(imV,ax=ax[1,0])
        ax[1,0].title.set_text('V')
        ax[1,0].set_xlim(left=Xpc[0,0],right=Xpc[0,-1])
        ax[1,0].set_ylim(bottom=Ypc[0,0],top=Ypc[-1,0])

        imP = ax[0,1].imshow(PdataResh, vmin=np.amin(Pdata), vmax=np.amax(Pdata), aspect='equal', origin='lower', extent= (Xpc[0,0], Xpc[0,-1]+dx, Ypc[0][0], Ypc[-1,0]+dy))
        cb = fig.colorbar(imP,ax=ax[0,1])
        ax[0,1].title.set_text('P')
        ax[0,1].set_xlim(left=Xpc[0,0],right=Xpc[0,-1])
        ax[0,1].set_ylim(bottom=Ypc[0,0],top=Ypc[-1,0])

    else:
        imU.set_data(UdataResh)
        imU.set_clim(vmin=np.amin(Udata),vmax=np.amax(Udata))

        imV.set_data(VdataResh)
        imV.set_clim(vmin=np.amin(Vdata),vmax=np.amax(Vdata))

        imP.set_data(PdataResh)
        imP.set_clim(vmin=np.amin(Pdata),vmax=np.amax(Pdata))

    plt.pause(0.001)

plt.show()

Ufile.close()
Vfile.close()
Pfile.close()
Rfile.close()
