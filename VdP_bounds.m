% -------------------------------------------------------------------------
% Ergodic optimization of the Van der Pol oscillator from data
%
% This script is used to replicate the results of Section 5.2 of Auxiliary
% Functions as Koopman Observables: Data-Driven Polynomial Optimization for
% Dynamical Systems by Jason J. Bramburger and Giovanni Fantuzzi.
%
% The goal of the script is to use data to approximate the Lie derivative 
% associated to the flow of the Van der Pol oscillator to bound long-time
% averages of the 'energy' g(x,x') = x^2 + x'^2. The energy is maximized 
% on the stable limit cycle, while it is minimized at the unstable
% steady-sate at (x,x') = (0,0).
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
% ---> See 
% maxPhi = max degree of phi dictionary of obserables
% maxPsi = max degree of psi dictionary of obserables
maxPhi = 10;
maxPsi = maxPhi + 2;

%% Generate synthetic data

% Integration parametes
dt = 0.001;
T = 10^(3); % integrate from t = 0 to t = T
t = 0:dt:T;
N = length(t);

% Initializations
x0 = [0.1; 0.2]; % initial condition with transients
% x0 = [1.325020770408186   1.573899618591706]; % starting on limit cycle
vdp = @(t,x) [x(2); -x(1) + 0.1*(1 - x(1)^2)*x(2)];

% Integrate system
[t, xdat] = ode45(vdp,t,x0,odeset('RelTol',1e-8,'AbsTol',1e-10));

% Plot trajectory
plot(xdat(:,1),xdat(:,2),'b','LineWidth',2)
xlabel('$x$','Interpreter','Latex')
ylabel('$\dot{x}$','Interpreter','Latex')
set(gca,'fontsize',16)

%% Compute empirical time average

obs = xdat(:,1).^2 + xdat(:,2).^2;
avg = trapz(t,obs)/t(end);

%% Create Psi_n and Phi_n matrices

% Phi matrix
pow = monpowers(2,maxPsi);
ell = size(pow,1); % number of nontrivial monomials in phi
Psi = zeros(ell,N-1);
for i = 1:ell
   zx = xdat(1:end-1,:).^pow(i,:);
   Psi(i,:) = prod(zx,2);
end

% Psi matrix
pow = monpowers(2,maxPhi);
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
x = sdpvar(2,1); % 2D state variable
z = monolist(x,maxPhi,0); % monomials that make up span(phi)
w = monolist(x,maxPsi,0); % monomials that make up span(psi)
c = sdpvar(length(z),1); % coeffs to build auxiliary function in span(phi)

% Lie approximation
% ---> Auxiliary function = c.'*z in span(phi)
L = (c.'*(K*w) - c.'*z)/dt;

%% Upper bound

sdpvar Bu 
cnstr1 = Bu - L - dot(x,x);
solvesos(sos(cnstr1),Bu,[],[c; Bu])

%% Lower bound

sdpvar Bl 
cnstr2 = L + dot(x,x) - Bl;
solvesos(sos(cnstr2),-Bl,[],[c; Bl])

%% Display results

fprintf('Numerically computed upper bound: %f \n',value(Bu))
fprintf('Emperical average computed from data: %f \n',avg)
fprintf('Numerically computed lower bound: %f \n',value(Bl))
fprintf('Exact lower bound: 0.0 \n')





