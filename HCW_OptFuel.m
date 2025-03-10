function [t_minU,X_minU,lam0] = HCW_OptFuel(x0,xf,tf,m0,lam0_guess,mu,T,c,rho,aT)
%%  < File Description >
%    Author:        Ruthvik Bommena
%    File Name:     HCW_OptFuel.m
%    Compiler:      MATLAB R2022b
%    Date:          22 March, 2023
%    Affiliation:   Department of Aerospace Engineering, University of Illinois Urbana-Champaign.
%    Description:   Function to propagate HCW states and costates differential equations for an Optimal Fuel Solution
%    Inputs:        Initial and final boundary conditions, time of flight, mass, costate guesses, 
%                   gravitational parameter, thrust, exhaust velocity, continuation parameter, target semi-major axis.

% Numerical Solution
opts_ode = odeset('RelTol',1e-13,'AbsTol',1e-15);
options = optimoptions('fsolve','Display','iter','MaxFunEvals',1e3,'MaxIter',1e3,'TolFun',1e-14,'TolX',1e-12,'UseParallel',false);

[lam0,~] = fsolve(@cost_minU,lam0_guess,options,0,tf,[x0; m0],xf,mu,T,c,rho,opts_ode,aT);

% Optimal Fuel Solution (min-U)
[t_minU,X_minU] = ode89(@HCW_OptF,[0 tf],[x0; m0; lam0],opts_ode,mu,T,c,rho,aT);

end

% Dynamics
function Xdot = HCW_OptF(~,X,mu,T,c,rho,aT)

x = X(1:6);
m = X(7);
lambda_r = X(8:10);
lambda_v = X(11:13); 
lambda_m = X(14);

S = norm(lambda_v)*c/m + lambda_m - 1; % Switch Function
delta = 0.5*(1+tanh(S/rho)); % Throttle

% accelaration components
ux = -T/m*delta*lambda_v(1)/norm(lambda_v); 
uy = -T/m*delta*lambda_v(2)/norm(lambda_v);
uz = -T/m*delta*lambda_v(3)/norm(lambda_v);

% State Dynamics
n = sqrt(mu/aT^3); % mean motion

xDot = [x(4);x(5);x(6);...
        3*n^2*x(1) + 2*n*x(5) + ux;...
       -2*n*x(4) + uy;...
       -n^2*x(3) + uz];

mDot = -T/c*delta;
     
% Costate Dynamics
H = zeros(3); H(1,2)=2*n; H(2,1)=-2*n; 
G = zeros(3); G(1,1) = 3*n^2; G(3,3) = -n^2;

lambda_rdot = -G'*lambda_v;
lambda_vdot = -lambda_r - H'*lambda_v;
lambda_mdot = -delta*T*norm(lambda_v)/(m^2);
          
Xdot = [xDot; mDot; lambda_rdot; lambda_vdot; lambda_mdot];

end

% Cost Function
function err = cost_minU(lam0_guess,t0,tf,x0,xf,mu,T,c,rho,opts_ode,aT)

[~,X] = ode89(@HCW_OptF,[t0 tf],[x0; lam0_guess],opts_ode,mu,T,c,rho,aT);

err = ([X(end,1:6),X(end,14)] - [xf;0]');

end
