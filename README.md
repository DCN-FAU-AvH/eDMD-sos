# **Data-Driven Polynomial Optimization for Dynamical Systems**

This repository contains MATLAB scripts to reproduce the data and figures from [Auxiliary Functions as Koopman Observables: Data-Driven Polynomial Optimization for Dynamical Systems](https://arxiv.org/abs/2303.01483) by [Jason J. Bramburger](https://hybrid.concordia.ca/jbrambur/) and [Giovanni Fantuzzi](https://dcn.nat.fau.eu/giovanni-fantuzzi/) (2023).

## **Paper Abstract**
The control Lyapunov function (CLF) approach to nonlinear control design is well established. Moreover, when the plant is control affine and polynomial, sum-of-squares (SOS) optimization can be used to find a polynomial controller as a solution to a semidefinite program. This letter considers the use of data-driven methods to design a polynomial controller by leveraging Koopman operator theory, CLFs, and SOS optimization. First, Extended Dynamic Mode Decomposition (EDMD) is used to approximate the Lie derivative of a given CLF candidate with polynomial lifting functions. Then, the polynomial Koopman model of the Lie derivative is used to synthesize a polynomial controller via SOS optimization. The result is a flexible data-driven method that skips the intermediary process of system identification and can be applied widely to control problems. The proposed approach is used to successfully synthesize a controller to stabilize an inverted pendulum on a cart.

## **Required MATLAB Packages**
All scripts require YALMIP and MOSEK to run. Both packages can be download for free at 
- YALMIP: https://yalmip.github.io/download/
- MOSEK: https://www.mosek.com/downloads/
