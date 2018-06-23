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
A warning message might show up at the end of a run on ptah as follows:
```BASH
[ptah:21834] 3 more processes have sent help message help-btl-vader.txt / cma-permission-denied
[ptah:21834] Set MCA parameter "orte_base_help_aggregate" to 0 to see all help / error messages
```
Then, add some MCA options. For instance,
```BASH
mpiexec  --mca orte_base_help_aggregate 0 --mca mca_btl_vader_single_copy_mechanism cma -np 4 ./laplace_MPI
```

[1]: https://www.psc.edu/images/xsedetraining/BootCamp2018/Boot_Camp_Hybrid_Computing_Challenge_copy.pdf
[2]: https://www.psc.edu/hpc-workshop-series/summer-bootcamp-2018
