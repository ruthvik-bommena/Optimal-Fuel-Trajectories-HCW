function [t_minU,X_minU,lam0,rho] = Solve_HCWOptFuel(x0,xf,tf,m0,mu,T,c,aT)
%%  < File Description >
%    Author:        Ruthvik Bommena
%    File Name:     Solve_HCWOptFuel.m
%    Compiler:      MATLAB R2022b
%    Date:          22 March, 2023
%    Affiliation:   Department of Aerospace Engineering, University of Illinois Urbana-Champaign.
%    Description:   Function to solve for the intial costate guesses for a min-fuel transfer and sweep rho to achieve bang-bang thrust profile
%    Inputs:        Initial and final boundary conditions, time of flight,mass, gravitational parameter, thrust,
%                   exhaust velocity, semi-major axis.

% Parameters
rho_f = 1e-4; err = 1; errTol = 1e-10; iter = 1; iterMax = 20;

% not a bang-bang solution
% rho = 1;
% lam0_guess = 1e-02*rand(7,1);

% converged rho and lambda for bang-bang solution
rho = 6.103515625000000e-05;
lam0_guess = [-1.302468642755131e-04;9.305255954692908e-04;1.590429190990787e-233;-2.922543204203493;4.392369611376025;5.432830646497902e-229;1.861057706982298e-05];

%% iterate
while err > errTol && iter < iterMax
    [t_minU,X_minU,lam0] = HCW_OptFuel(x0,xf,tf,m0,lam0_guess,mu,T,c,rho,aT);
    err = norm([X_minU(end,1:6),X_minU(end,14)] - [xf;0]');
    lam0_guess = lam0;
    iter = iter+1; 
end

% Sweep rho to get bang-bang thrust profile
while rho > rho_f
    rho = rho*0.5;
    [t_minU,X_minU,lam0] = HCW_OptFuel(x0,xf,tf,m0,lam0,mu,T,c,rho,aT);
end

end 