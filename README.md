# **Data-Driven Polynomial Optimization for Dynamical Systems**

This repository contains MATLAB scripts to reproduce the data and figures from [Auxiliary Functions as Koopman Observables: Data-Driven Polynomial Optimization for Dynamical Systems](https://arxiv.org/abs/2303.01483) by [Jason J. Bramburger](https://hybrid.concordia.ca/jbrambur/) and [Giovanni Fantuzzi](https://dcn.nat.fau.eu/giovanni-fantuzzi/) (2023).

## **Paper Abstract**
We present a flexible data-driven method for dynamical system analysis that does not require explicit model discovery. The method is rooted in well-established techniques for approximating the Koopman operator from data and is implemented as a semidefinite program that can be solved numerically. Furthermore, the method is agnostic of whether data is generated through a deterministic or stochastic process, so its implementation requires no prior adjustments by the user to accommodate these different scenarios. Rigorous convergence results justify the applicability of the method, while also extending and uniting similar results from across the literature. Examples on discovering Lyapunov functions, performing ergodic optimization, and bounding extrema over attractors for both deterministic and stochastic dynamics exemplify these convergence results and demonstrate the performance of the method. 

## **Required MATLAB Packages**
All scripts require YALMIP and MOSEK to run. Both packages can be download for free at 
- YALMIP: https://yalmip.github.io/download/
- MOSEK: https://www.mosek.com/downloads/
