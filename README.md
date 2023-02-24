# **Data-Driven Polynomial Optimization for Dynamical Systems**

This repository contains MATLAB scripts to reproduce the data and figures from [Auxiliary Functions as Koopman Observables: Data-Driven Polynomial Optimization for Dynamical Systems](TBD) by Jason J. Bramburger and Giovanni Fantuzzi (2023).

## **Paper Abstract**
We present a flexible data-driven method for dynamical system analysis that does not require explicit model discovery. The method is rooted in well-established tech- niques for approximating the Koopman operator from data and is implemented as a semidefinite program that can be solved numerically. The method is agnostic of whether data is generated through a deterministic or stochastic process, so its imple- mentation requires no prior adjustments by the user to accommodate these different scenarios. Rigorous convergence results justify the applicability of the method, while also extending and uniting similar results from across the literature. Examples on dis- covering Lyapunov functions and on performing ergodic optimization for both deter- ministic and stochastic dynamics exemplify these convergence results and demonstrate the performance of the method.

## **Required Packages**
All scripts require YALMIP and MOSEK to run. Both packages can be download for free at 
- YALMIP: https://yalmip.github.io/download/
- MOSEK: https://www.mosek.com/downloads/

The script logistic_bounds.m further uses the ChebFun package to improve numerical performance. It can be download for free at: https://www.chebfun.org/download/ 

