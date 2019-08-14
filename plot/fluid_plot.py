#!/usr/bin/env python3.7
import numpy as np
import seaborn
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import pandas as pd

# Read specific parameters
def read_input(filename, var):
    inData = pd.read_csv(inputfile,delimiter= '=', names = ["variable", "value"])
    return_val = np.empty(np.shape(var))
    for i, var_val in enumerate(var):
        val = inData[inData.variable.str.strip() == var_val].values[0,1]
        return_val[i] = val
    return return_val[:]

# Create a mesh for the output of the simulation
def make_pcolor_mesh(xl,yl,im,jm):
    dx = xl/im
    dy = yl/jm
    x = np.arange(0,xl,dx)
    y = np.arange(0,yl,dy)
    # This is usual so that we can see the last data
    x = np.append(x,x[-1]+dx)
    y = np.append(y,y[-1]+dy)
    # Add a cell before and after (for the boundaries added in the simulation)
    x = np.insert(x,0,x[0]-dx)
    x = np.append(x,x[-1]+dx)
    y = np.insert(y,0,y[0]-dy)
    y = np.append(y,y[-1]+dy)
    # Create mesh
    X, Y = np.meshgrid(x,y)
    return(X,Y)

# Convert the imshow mesh into the positions for quiver (i.e. center of cells)
def pcolor_mesh_to_quiver(X,Y,dx,dy):
    Xq = X[:-1,:-1]+dx/2
    Yq = Y[:-1,:-1]+dy/2
    return(Xq,Yq)

# Interpolate vel fields to center of cells
def interpolate_data(U):
    U1 = U[:,1:,1:]
    U2 = U[:,:-1,:-1]
    Uout = (U1+U2)/2
    return Uout

# Read all the output data and reshape them so that dim = [time,jmax+2,imax+2]
def read_output_lines(filename,l,m):
    file = open(filename, "r")
    lines = file.read().strip().split("\n")
    out = np.empty((len(lines),m,l))
    for i,line in enumerate(lines):
        data = np.fromstring(line, dtype=float, sep=' ')
        out[i,:,:] = np.reshape(data,(m,l))
    file.close()
    return out

def time_evolution_imshow(X,Y,U,V,P):
    # Append useless last lines for plot
    newshape = [x+1 if i > 0 else x for i,x in enumerate(U.shape)]
    Udata = np.zeros(newshape)
    Udata[:,:-1,:-1] = U
    Vdata = np.zeros(newshape)
    Vdata[:,:-1,:-1] = V
    Pdata = np.zeros(newshape)
    Pdata[:,:-1,:-1] = P

    fig, ax = plt.subplots(2, 2)
    for i in range(newshape[0]):
        if (i==0):
            imU = ax[0,0].imshow(Udata[i,:,:], vmin=np.amin(Udata[i,:,:]), vmax=np.amax(Udata[i,:,:]),
                    aspect='equal', origin='lower', extent= (X[0,0], X[0,-1]+(X[0,1]-X[0,0]), Y[0][0], Y[-1,0]+(X[1][0])-X[0][0]))
            cb = fig.colorbar(imU,ax=ax[0,0])
            ax[0,0].title.set_text('U')
            ax[0,0].set_xlim(left=X[0,0],right=X[0,-1])
            ax[0,0].set_ylim(bottom=Y[0,0],top=Y[-1,0])

            imV = ax[1,0].imshow(Vdata[i,:,:], vmin=np.amin(Vdata[i,:,:]), vmax=np.amax(Vdata[i,:,:]),
                    aspect='equal', origin='lower', extent= (X[0,0], X[0,-1]+(X[0,1]-X[0,0]), Y[0][0], Y[-1,0]+(X[1][0])-X[0][0]))
            cb = fig.colorbar(imV,ax=ax[1,0])
            ax[1,0].title.set_text('V')
            ax[1,0].set_xlim(left=X[0,0],right=X[0,-1])
            ax[1,0].set_ylim(bottom=Y[0,0],top=Y[-1,0])

            imP = ax[0,1].imshow(Pdata[i,:,:], vmin=np.amin(Pdata[i,:,:]), vmax=np.amax(Pdata[i,:,:]),
                    aspect='equal', origin='lower', extent= (X[0,0], X[0,-1]+(X[0,1]-X[0,0]), Y[0][0], Y[-1,0]+(X[1][0])-X[0][0]))
            cb = fig.colorbar(imP,ax=ax[0,1])
            ax[0,1].title.set_text('P')
            ax[0,1].set_xlim(left=X[0,0],right=X[0,-1])
            ax[0,1].set_ylim(bottom=Y[0,0],top=Y[-1,0])
        else:
            imU.set_data(Udata[i,:,:])
            imU.set_clim(vmin=-0.02,vmax=0.02)#np.amin(Udata[i,:,:]),vmax=np.amax(Udata[i,:,:]))

            imV.set_data(Vdata[i,:,:])
            imV.set_clim(vmin=-0.02,vmax=0.02)#np.amin(Vdata[i,:,:]),vmax=np.amax(Udata[i,:,:]))

            imP.set_data(Pdata[i,:,:])
            imP.set_clim(vmin=np.amin(Pdata[i,:,:]),vmax=np.amax(Pdata[i,:,:]))
        plt.pause(0.01)
    plt.show()

def time_evolution_quiver(X,Y,U,V):
    res = 2
    # Append useless last lines for plot
    newshape = [x+1 if i > 0 else x for i,x in enumerate(U.shape)]
    Udata = np.zeros(newshape)
    Udata[:,:-1,:-1] = U
    Vdata = np.zeros(newshape)
    Vdata[:,:-1,:-1] = V
    UN = np.sqrt( Udata**2 + Vdata**2 )
    Udata = np.divide(Udata,UN)
    Vdata = np.divide(Vdata,UN)

    fig, ax = plt.subplots(1,1)
    for i in range(newshape[0]):
        if (i==0):
            q = ax.quiver(X[::res,::res],Y[::res,::res],Udata[i,::res,::res], Vdata[i,::res,::res], UN[i,::res,::res])
            cb = fig.colorbar(q,ax=ax)
            ax.title.set_text('Velocity field')
            #ax.set_xlim(left=X[0,0],right=X[0,-1])
            ax.set_xlim(left=5,right=15)
            ax.set_ylim(bottom=Y[0,0],top=Y[-1,0])
            ax.set_aspect('equal', adjustable='box')
        else:
            q.set_UVC(Udata[i,::res,::res],Vdata[i,::res,::res], UN[i,::res,::res])
            q.set_clim(vmin=np.amin(UN[i,::res,::res]),vmax=np.amax(UN[i,::res,::res]))

        plt.pause(0.01)
    plt.show()

#------------------------------

#--------------------------------------------------
#       Start of the program
#--------------------------------------------------
inputfile = "../IO/input"
outputfile = "../IO/output"

var = ["xlength","ylength","imax","jmax"];

[xlength,ylength,imax,jmax] = read_input(inputfile,var)
imax = int(imax)
jmax = int(jmax)
dx = xlength/imax
dy = ylength/jmax

(Xpc,Ypc) = make_pcolor_mesh(xlength,ylength,imax,jmax)

(Xq,Yq) = pcolor_mesh_to_quiver(Xpc,Ypc,dx,dy)

# # U Defined on the middle of the right edge of each cells
# XU = Xpc[:-1,:-1] + dx
# YU = Ypc[:-1,:-1] + dy/2
# # V Defined on the middle of the top edge of each cells
# XV = Xpc[:-1,:-1] + dx/2
# YV = Ypc[:-1,:-1] + dy
# # V Defined in the middle of each cells
# XP = Xpc[:-1,:-1] + dx/2
# YP = Ypc[:-1,:-1] + dy/2

Udata = read_output_lines(outputfile+'U',imax+2,jmax+2)
Vdata = read_output_lines(outputfile+'V',imax+2,jmax+2)
Pdata = read_output_lines(outputfile+'P',imax+2,jmax+2)

Ucenter = interpolate_data(Udata)
Vcenter = interpolate_data(Vdata)

#time_evolution_imshow(Xpc,Ypc,Udata,Vdata,Pdata)
time_evolution_quiver(Xq,Yq,Ucenter,Vcenter)

#   #------------------------------
#   #   Make gif
#   import sys
#   from matplotlib.animation import FuncAnimation
#   X = Xq
#   Y = Yq
#   U = Ucenter
#   V = Vcenter
#   newshape = [x+1 if i > 0 else x for i,x in enumerate(U.shape)]
#   Udata = np.zeros(newshape)
#   Udata[:,:-1,:-1] = U
#   Vdata = np.zeros(newshape)
#   Vdata[:,:-1,:-1] = V
#   UN = np.sqrt( Udata**2 + Vdata**2 )
#
#   fig, ax = plt.subplots()
#   fig.set_tight_layout(True)
#
#   print('fig size: {0} DPI, size in inches {1}'.format(
#           fig.get_dpi(), fig.get_size_inches()))
#   res = 2
#   q = ax.quiver(X[::res,::res],Y[::res,::res],Udata[0,::res,::res], Vdata[0,::res,::res], UN[0,::res,::res])
#   cb = fig.colorbar(q,ax=ax)
#   ax.title.set_text('Velocity field')
#   ax.set_xlim(left=X[0,0],right=X[0,-1])
#   ax.set_ylim(bottom=Y[0,0],top=Y[-1,0])
#   ax.set_aspect('equal', adjustable='box')
#
#
#   def update(i):
#       q.set_UVC(Udata[i,::res,::res],Vdata[i,::res,::res], UN[i,::res,::res])
#       q.set_clim(vmin=np.amin(UN[i,::res,::res]),vmax=np.amax(UN[i,::res,::res]))
#       return q
#   anim = FuncAnimation(fig, update, frames=np.arange(0, 40), interval=200)
#   anim.save('quiver.gif', dpi=200, writer='imagemagick')
#   #------------------------------
