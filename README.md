# **Data-Driven Polynomial Optimization for Dynamical Systems**

This repository contains MATLAB scripts to reproduce the data and figures from [Auxiliary Functions as Koopman Observables: Data-Driven Polynomial Optimization for Dynamical Systems](https://link.springer.com/article/10.1007/s00332-023-09990-2) by [Jason J. Bramburger](https://hybrid.concordia.ca/jbrambur/) and [Giovanni Fantuzzi](https://dcn.nat.fau.eu/giovanni-fantuzzi/) (2023).

## **Paper Abstract**
We present a flexible data-driven method for dynamical system analysis that does not require explicit model discovery. The method is rooted in well-established techniques for approximating the Koopman operator from data and is implemented as a semidefinite program that can be solved numerically. Furthermore, the method is agnostic of whether data is generated through a deterministic or stochastic process, so its implementation requires no prior adjustments by the user to accommodate these different scenarios. Rigorous convergence results justify the applicability of the method, while also extending and uniting similar results from across the literature. Examples on discovering Lyapunov functions, performing ergodic optimization, and bounding extrema over attractors for both deterministic and stochastic dynamics exemplify these convergence results and demonstrate the performance of the method. 

## **Required MATLAB Packages**
All scripts require YALMIP and MOSEK to run. Both packages can be download for free at 
- YALMIP: https://yalmip.github.io/download/
- MOSEK: https://www.mosek.com/downloads/

## **Repository Contents**
This repository contains MATLAB script to reproduce the results for the examples in Section 6 of the paper. Precisely, the scripts perform the following tasks:
- `lyapunov_function.m` discovers a Lyapunov function from data and corresponds to the example in Section 6.1.
- `VdP_bounds.m` provides upper and lower bounds on the long-time average of the energy for the Van der Pol oscillator. This script accompanies the example in Section 6.2 and, in particular, is used to generate the data in Table 1.
- `logistic_bounds.m` provides upper and lower bounds on the long-time average of the state-space observable, x, for the random logistic map. Chebyshev basis functions are used to improve numerical conditioning and therefore one must include the ChebFun package (see above for download link). The MATLAB functions `chebsdp_1d.m`, `chebsdp_1d_locball.m`, and `edmd_with_thresholding.m` are all called by `logistic_bounds.m` to perform specific tasks, primarily related to the Chebyshev basis functions. This script accompanies the example in Section 6.3 and, in particular, is used to generate the data in Table 2.
- `thomas_bounds.m` provides pointwise bounds on [Thomas' cyclically symmetric attractor](https://en.wikipedia.org/wiki/Thomas%27_cyclically_symmetric_attractor) using a data-driven implementation of the methods of [Goluskin (2020)](https://iopscience.iop.org/article/10.1088/1361-6544/ab8f7b). This script corresponds to the example in Section 6.4 and is used to produce the data in Table 3 of the manuscript.  
