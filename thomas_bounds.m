% -------------------------------------------------------------------------
% Pointwise attractor bounds on Thomas' cyclically symmetric attractor
%
% This script is used to replicate the results of Section 6.4 of Auxiliary
% Functions as Koopman Observables: Data-Driven Polynomial Optimization for
% Dynamical Systems by Jason J. Bramburger and Giovanni Fantuzzi.
%
% The goal of the script is to use implement a data-driven formulation of
% Goluskin's pointwise attractor bounding method (Nonlinearity, 2018) on
% Thomas' cyclically symmetric system:
%
%          x_i' = sin(x_{i+1}) - 0.2*x_i,      i = 1,2,3,
%
% where x_4 = x_1 for the presentation. 
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
maxPhi = 4;
maxPsi = 4;

%% Generate synthetic data

% Integration parametes
dt = 0.001;
T = 1e2; % integrate from t = 0 to t = T
t = 0:dt:T;
N = length(t);

% System parameter (as b --> 0 system becomes more chaotic)
b = 0.2;

% Integrate system
x0 = [0.1; 0.2;0.3]; 
thomas = @(t,x) [sin(x(2)) - b*x(1); sin(x(3)) - b*x(2); sin(x(1)) - b*x(3)]; 
[t, xdat] = ode45(thomas,t,x0);

% Scale data into hypercube [-1,1]^3
scl = 5;
xdat = xdat/scl;

% Plot trajectory
plot3(xdat(:,1),xdat(:,2),xdat(:,3),'b','LineWidth',2)
xlabel('$x_1$','Interpreter','Latex')
ylabel('$x_2$','Interpreter','Latex')
zlabel('$x_3$','Interpreter','Latex')
set(gca,'fontsize',16)
grid on

%% Create Psi_n and Phi_n matrices

% Phi matrix
pow = monpowers(3,maxPsi);
ell = size(pow,1); % number of nontrivial monomials in phi
Psi = zeros(ell,N-1);
for i = 1:ell
   zx = xdat(1:end-1,:).^pow(i,:);
   Psi(i,:) = prod(zx,2);
end

% Psi matrix
pow = monpowers(3,maxPhi);
m = size(pow,1); % number of nontrivial monomials in phi
Phi = zeros(m,N-1);
for i = 1:m
   zy = xdat(2:end,:).^pow(i,:);
   Phi(i,:) = prod(zy,2)';
end

%% Create Koopman and Lie derivative matrix approximations

% Koopman approximation
K = Phi*pinv(Psi);

% Symbolic variables
x = sdpvar(3,1); % 3D state variable
z = monolist(x,maxPhi,0); % monomials that make up span(phi)
w = monolist(x,maxPsi,0); % monomials that make up span(psi)
c = sdpvar(length(z),1); % coeffs to build auxiliary function in span(phi)

% Lie approximation
% ---> Auxiliary function = c.'*z in span(phi)
thresh = 1e-4; % threshold out small values for interpretability
Lie = (K - eye(size(K)))/dt;
Lie(abs(Lie) < thresh) = 0;
L = c.'*Lie*w;
%L = (c.'*(K*w) - c.'*z)/dt;

%% Attractor bounds

% S procedure
[s1,sc1] = polynomial(x,ceil(maxPhi/2)*2 - 2);
[s2,sc2] = polynomial(x,ceil(maxPhi/2)*2 - 2);

% Optimization procedure parameters
sdpvar C % optimization objective
lam = 1;

% x1 bound
cnstr = [sos(c.'*z - x(1) - (1 - x(1)^2 - x(2)^2 - x(3)^2)*s1); sos(C - c.'*z - lam*L - (1 - x(1)^2 - x(2)^2 - x(3)^2)*s2), sos(s1), sos(s2)];
solvesos(cnstr,C,[],[c;C;sc1;sc2])
x1bnd = (scl)*value(C);

% x2 bound
cnstr = [sos(c.'*z - x(2) - (1 - x(1)^2 - x(2)^2 - x(3)^2)*s1); sos(C - c.'*z - lam*L - (1 - x(1)^2 - x(2)^2 - x(3)^2)*s2), sos(s1), sos(s2)];
solvesos(cnstr,C,[],[c;C;sc1;sc2])
x2bnd = (scl)*value(C);


% x3 bound
cnstr = [sos(c.'*z - x(3) - (1 - x(1)^2 - x(2)^2 - x(3)^2)*s1); sos(C - c.'*z - lam*L - (1 - x(1)^2 - x(2)^2 - x(3)^2)*s2), sos(s1), sos(s2)];
solvesos(cnstr,C,[],[c;C;sc1;sc2])
x3bnd = (scl)*value(C);

% Print results
fprintf('Numerically computed upper bound on x_1: %f \n',x1bnd)
fprintf('Numerically computed upper bound on x_2: %f \n',x2bnd)
fprintf('Numerically computed upper bound on x_3: %f \n',x3bnd)




