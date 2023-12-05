function [t_minU,X_minU,lam0,rho] = Solve_HCWOptFuel(x0,xf,tf,m0,mu,T,c,aT)
%%  < File Description >
%    Author:        Ruthvik Bommena
%    File Name:     Solve_HCWOptFuel.m
%    Compiler:      MATLAB R2022b
%    Date:          22 March, 2023
%    Affiliation:   Department of Aerospace Engineering, University of Illinois Urbana-Champaign.
%    Description:   Function to solve minU transfer and sweep rho to achieve bang-bang thrust profile
%    Inputs:        Initial and final boundary conditions, time of flight,mass, gravitational parameter, thrust,
%                   exhaust velocity, semi-major axis.

% Parameters
rho_f = 1e-4; err = 1; errTol = 1e-10; iter = 1; iterMax = 20;

% rho = 1;
% lam0_guess = 1e-02*rand(7,1);

% rho and lambda guesses
rho = 9.404610869860067e-05;
lam0_guess = [-0.017172802398584;5.970717937722085e-04;1.952653584036095e-82;-1.658252884212246;-8.013680679648221;7.066155254214553e-80;1.196168943379431e-05];

% Solve Optimal-Fuel transfer
while err > errTol && iter < iterMax
    [t_minU,X_minU,lam0] = HCW_OptFuel(x0,xf,tf,m0,lam0_guess,mu,T,c,rho,aT);
    err = norm([X_minU(end,1:6),X_minU(end,14)] - [xf;0]');
    lam0_guess = lam0;
    iter = iter+1; 
end

% Sweep rho to get bang-bang thrust profile
while rho > rho_f
    rho = rho*0.9;
    [t_minU,X_minU,lam0] = HCW_OptFuel(x0,xf,tf,m0,lam0,mu,T,c,rho,aT);
end

end 