# laplace_chlng

This code computes a static equilibrium temperature distribution in a two-dimensional domain by approximately solving the Laplace equation with boundary conditions using the finite difference method.


## How to build
```BASH
mpicc ./laplace_mpi_sol.c -o ./laplace_mpi_sol
```

## How to run
```BASH
mpiexec -np 4 ./laplace_mpi_sol
```
