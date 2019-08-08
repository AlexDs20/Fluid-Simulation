# An incompressible fluid simulation using finite difference method

This simulation is implemented in C++ and follows the description given in Numerical simulation in Fluid Dynamics A Practical Introduction by M. Griebel, T. Dornseifer and T. Neunhoeffer.

### What is implemented
- [x] 2D simulation in a rectangular grid
- [x] Incompressible fluid
- [x] 5 different boundary conditions at the edges of the grid: no-slip, free-slip,
- [x] outflow, inflow, periodic boundaries
- [x] A python script to visualise the output

### What is *not* implemented

- [ ] Interpolation of the output data to make a correct quiver plot
- [ ] Obstacles inside the domain
- [ ] Obstacles moving due to the flow
- [ ] 3D simulation
- [ ] Parallelized computation

### What is required

All of this was developed on a linux machine and may not work directly on
Windows or Mac.

* mpich package (openMPI, g++) for compiling and running (even though parallel computation is not implemented yet).
I use version 3.3.1-1 of mpich found on the AUR.

* The libraries numpy, matplotlib and pandas in python to visualise the output. I use python 3.7

## A typical run

A typical run would consist of 3:
1) Set the inputs for the desired simulation
2) Run the simulation
3) Look at the outputs

### The inputs

The input file is located at IO/input:

Inputs     | what it does
-----------|--------------
*xlength*  |   length of the system along x
*ylength*  |   length of the system along y
*imax*     |   number of cells in x
*jmax*     |   number of cells in y
*dt*       |   minimum desired time step
*t_end*    |   simulation time
*tau*      |   security factor (0<tau<=1) for the time step
*itermax*  |   maximum number of iteration to get the pressure when using the SOR algo.
*eps*      |   maximum residual in the SOR algo.
*omega*    |   coeff in the SOR algo
*Re*       |   Reynolds number
*gx*       |   external force along x
*gy*       |   external force along y
*UI*       |   Initial velocity along x (given in the whole domain)
*VI*       |   Initial velocity along y (given in the whole domain)
*PI*       |   Initial pressure (given in the whole domain)
*rho*      |   density
*dynvis*   |   dynamical viscosity
*wW*       |   boundary condition on the left of the box (x=0)
*wE*       |   boundary condition on the right of the box (x=xlength)
*wN*       |   boundary condition on the top of the box (y=ylength)
*wS*       |   boundary condition on the bottom of the box (y=0)


The values of *rho* and *dynvis* are not currently used.
I hesitated between requesting physical inputs values or dimensionless inputs.

Currently the inputs should be dimensionless and the parameters *rho* and *dynvis* are hence useless.


The *wW*, *wE*, *wS*, *wN* can take 5 values: 1, 2, 3 ,4 or 5.
A value of:
* 1: indicates the no slip condition i.e. all the velocities vanish at the boundary
* 2: indicates the free slip condition i.e. normal velocity is 0 but the tangential is "free" ($ \frac{\partial \phi_t}{\partial n} = 0 $).
* 3: indicates outflow i.e. the fluid simply gets out of the domain.
* 4: indicates inflow i.e. fluid enters the domain. (There is no current input for the velocity of the flow, it is hard coded in the boundary.cpp file)
* 5: indicates periodic boundaries. For this one to work, the pair wW and wE or the pair wS and wN should be given the value 5.


If you would like more specific boundary conditions, such as a narrow inflow,
this needs to be implemented in the setSpecBoundaries function in the
boundary.cpp file.
Doing this will require recompiling the code.

### To run

If you edited the code, run:
```
make
```

You can now run the code:
```
make run
```


### To visualise the output

Go the the plot directory and run:
```
./fluid_plot.py
```

Note that the output of the simulation of the velocities are at the border of the cells,
to make a quiver plot starting from the center of the cells, the values of the velocity fields must be evaluated there first.
