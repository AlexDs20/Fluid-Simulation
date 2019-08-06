# An incompressible fluid simulation using finite difference method

This simulation is implemented in C++ and follows the description given in Numerical simulation in Fluid Dynamics A Practical Introduction by M. Griebel, T. Dornseifer and T. Neunhoeffer.

Currently, it is only in 2D and only the boundary conditions at the edges of the box have been implemented.
Even though it requires mpich for compiling, nothing is parallelized yet.

## To run
```
make clean
make 
make run
```

The input for the simulation are given in IO/input:

* xlength      length of the system along x
* ylength      length of the system along y
* imax         number of cells in x
* jmax         number of cells in y
* dt           minimum desired time step
* t_end        simulation time
* tau          security factor (0<tau<=1) for the time step
* itermax      maximum number of iteration to get the pressure when using the SOR algo.
* eps          maximum residual in the SOR algo.
* omega        coeff in the SOR algo
* Re           Reynolds number
* gx           external force along x
* gy           external force along y                        
* UI           Initial velocity along x (given in the whole domain)
* VI           Initial velocity along y (given in the whole domain)
* PI           Initial pressure (given in the whole domain)
* rho          density
* dynvis       dynamical viscosity
* wW           boundary condition on the left of the box (x=0)
* wE           boundary condition on the right of the box (x=xlength)
* wN           boundary condition on the top of the box (y=ylength)
* wS           boundary condition on the bottom of the box (y=0)




The values of *rho* and *dynvis* are not currently used.

I hesitated between requesting physical inputs values or dimensionless inputs.

Currently the inputs should be dimensionless and these parameters are hence useless.



The *wW*, *wE*, *wS*, *wN* can take 5 values: 1, 2, 3 ,4 or 5.
A value of:
* 1: indicates the no slip condition i.e. all the velocities vanish at the boundary
* 2: indicates the free slip condition i.e. normal velocity is 0 but the tangential is free.
* 3: indicates outflow i.e. the fluid simply gets out of the domain.
* 4: indicates inflow i.e. fluid enters the domain. (There is no current input for the velocity of the flow, it is hard coded in the boundary.cpp file)
* 5: indicates periodic boundaries. For this one to work, the pair wW and wE or the pair wS and wN should be given the value 5.



## To visualise

Go the the plot directory and run:
```
./fluid_plot.py
```

Note that the output of the simulation of the velocities are at the border of the cells,
to make a quiver plot starting from the center of the cells, the values of the velocity fields must be evaluated there first.

## On the to do list:

1) Implementing obstacle within the domain so that the user can simulate whatever he wants.
2) Paralellizing the code.
3) Making the possibility for the obstacle to be moved by the flows.
4) Giving the possibility to do 3D simulations
5) Making the quiver plot to visualise the velocity field.
