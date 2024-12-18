# Finite element simulation of the deterministic PDEs describing coarse grained biased QSAPs

## Installation

The code has been tested with fenicsX v0.7.3 alongside the dolfinX multi-point constraints package v0.7.2. To install these and other required dependencies in a conda environment called `AMIPS-env` enter the following commands

```
conda create --name AMIPS-env python=3.11
conda activate AMIPS-env
conda install -c conda-forge fenics-dolfinx=0.7.3 dolfinx_mpc=0.7.2 matplotlib
```

## Running the code

To run the code over (for example) 8 threads with $\beta = 1/2$, $\text{Pe}=7$, and $\rho_0=0.7$ run

`mpirun -n 8  python3 ./run_amips.py 0.5 7.0 0.7 -m exampleDirectorySuffix`

Other parameters and initial conditions can be set in the `./run_amips.py` file.

If the selected model is the full system in $\{\rho,\psi\}$, this will create a directory `./data/beta_0.5_Pe_7.0_p_0.7_full_exampleDirectorySuffix`

If the selected model is the implict nsAMB model in $\phi$, this will create a directory `./data/beta_0.5_Pe_7.0_p_0.7_NSAMB_exampleDirectorySuffix`