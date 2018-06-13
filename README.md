# laplace_chlng

This code computes a static equilibrium temperature distribution in a square domain by approximately solving the Laplace equation with boundary conditions using the finite difference method.

Optimally accelerating this code using any combination of OpenMP, OpenACC and MPI is the goal of [the 5th Annual XSEDE Boot Camp Hybrid Challenge][1]. This challenge is a part of the activities of [XSEDE HPC Workshop: Summer Boot Camp][2] held at the Pittsburgh Supercomputing Center, June 4-7, 2018.


## How to build
```BASH
mpicc ./laplace_mpi_sol.c -o ./laplace_mpi_sol
```

## How to run
```BASH
mpiexec -np 4 ./laplace_mpi_sol
```

[1]: https://www.psc.edu/images/xsedetraining/BootCamp2018/Boot_Camp_Hybrid_Computing_Challenge_copy.pdf
[2]: https://www.psc.edu/hpc-workshop-series/summer-bootcamp-2018
