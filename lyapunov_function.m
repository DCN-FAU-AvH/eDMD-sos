% -------------------------------------------------------------------------
% Discovering Lyapunov functions from data
%
% This script is used to replicate the results of Section 6.1 of Auxiliary
% Functions as Koopman Observables: Data-Driven Polynomial Optimization for
% Dynamical Systems by Jason J. Bramburger and Giovanni Fantuzzi.
%
% The goal of the script is to use data to approximate the Lie derivative 
% associated to the flow of the dynamical system and then use the
% approximate Lie derivative to discover a Lyapunov function that provides
% the global stability of the system associated to the data. At the end of
% the script we further confirm that the discovered function is indeed a 
% Lyapunov for the system using the exact Lie derivative from the known 
% system dynamics.
%
% Packages required: YALMIP and MOSEK
%
% Written by J. Bramburger and G. Fantuzzi.
%
% -------------------------------------------------------------------------

% Clean workspace
clear; close all; clc
yalmip clear
format long

%% Method Parameters 
% maxPhi = max degree of phi dictionary of obserables
% maxPsi = max degree of psi dictionary of obserables
% epsilon = hyperparameter specific to Lyapunov function for sharp bounds
maxPhi = 4;
maxPsi = 8;
epsilon = 1;

%% Generate synthetic data

% Number of data points
N = 1e4; 

% Random points in [-2,2]x[-2,2]
xdat = 2*rand(N,2) - 1;
ydat = zeros(N,2);

% Images of random points under map dynamics
ydat(:,1) = 0.3*xdat(:,1);
ydat(:,2) = -xdat(:,1) + 0.5*xdat(:,2) + (7/18)*xdat(:,1).^2;

%% Create Psi_n and Phi_n matrices

% Psi matrix
pow = monpowers(2,maxPsi);
ell = size(pow,1); % number of nontrivial monomials in psi
Psi = zeros(ell,N);
for i = 1:ell
   zx = xdat.^pow(i,:);
   Psi(i,:) = prod(zx,2);
end

% Phi matrix
pow = monpowers(2,maxPhi);
m = size(pow,1); % number of nontrivial monomials in phi
Phi = zeros(m,N);
for i = 1:m
   zy = ydat.^pow(i,:);
   Phi(i,:) = prod(zy,2)';
end

%% Create Koopman and Lie derivative matrix approximations

% Koopman approximation
K = Phi*pinv(Psi);

% Symbolic variables
x = sdpvar(2,1); % 2D state variable
z = monolist(x,maxPhi,0); % monomials that make up span(phi)
w = monolist(x,maxPsi,0); % monomials that make up span(psi)
c = sdpvar(length(z),1); % coeffs to build Lyapunov function in span(phi)

% Lie approximation
% ---> Lyapunov function = c.'*z in span(phi)
L = c.'*(K*w) - c.'*z;

%% Identify Lyapunov function with SOS programming

% Inequality constraints posed relaxed to SOS conditions
cons = [sos(c.'*z - epsilon*dot(x,x)); sos(-L - epsilon*dot(x,x))];

% Set solver to MOSEK
opts = sdpsettings;
opts.solver = 'mosek';

% Objective function: minimize l^1 norm of coefficients
OBJ = sum(abs(c));

% Solve SOS problem
solvesos(cons,OBJ,opts,c)

% Print coefficients after removing those smaller than 10^-4 
c = clean(value(c), 1e-4)

%% Check identified Lyapunov function is a real Lyapunov function
% ---> This time we maximize epsilon and check if it is positive

% Initializations
yalmip clear
x = sdpvar(2,1); % 2D state variable
z = monolist(x,maxPhi,0); % monomials that make up span(phi)
sdpvar epsilon % epsilon is now a variable to be maximized

% Discovered Lyapunov function
v = c.'*z;

% True right-hand-side of map
f = [0.3*x(1); -x(1) + 0.5*x(2) + (7/18)*x(1)^2];

% True Lie derivative
Lie_exact = replace(v,x,f);

% Find Lyapunov function
solvesos([sos(v - epsilon*dot(x,x)); sos(- (Lie_exact - v) - epsilon*dot(x,x))],-epsilon)
fprintf('Maximized value of epsilon is: %f \n',value(epsilon))

if value(epsilon) > 0
   fprintf('We have discovered a true Lyapunov function from the data! \n') 
else
    fprintf('Something went wrong - we cannot verify that this is a Lyapunov function. \n')
end

 





