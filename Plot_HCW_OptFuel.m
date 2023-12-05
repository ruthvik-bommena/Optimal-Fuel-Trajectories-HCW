function Plot_HCW_OptFuel
%%  < File Description >
%    Author:        Ruthvik Bommena
%    File Name:     Plot_HCW_OptFuel.m
%    Compiler:      MATLAB R2022b
%    Date:          22 March, 2023
%    Affiliation:   Department of Aerospace Engineering, University of Illinois Urbana-Champaign.
%    Description:   Function to plot HCW optimal fuel transfer for chaser to rendezvous with the target.

close all; clear; clc;

% Parameters
aT = 6798.1; % km
mu = 398600; % km^3/s^2
m0 = 50; % mass = 50 kg
T =  0.0015/1000; % MiXi 1.5mN
Isp = 4190; % s
g0 = 9.8/1000; % km/s^2
c = Isp*g0/1000; % km/s

% Boundary conditions
x0 = [0 1 0 0 0 0]'; % km-s
xf = [0 0 0 0 0 0]'; % km-s

% Transfer time
tf = 120; % min
tf = tf*60; % sec

% OptFuel transfer
[t_minU,X_minU,~,rho] = Solve_HCWOptFuel(x0,xf,tf,m0,mu,T,c,aT);

% Plots
plots(t_minU,X_minU,rho,m0)

end


%% Function Plots
function plots(t_minU,X_minU,rho,m0)

% Get relative distance and velocity
position_err = zeros(length(X_minU),1);
velocity_err = zeros(length(X_minU),1);
for ii = 1:length(X_minU)
    position_err(ii) = norm(X_minU(ii,1:3))*1000;
    velocity_err(ii) = norm(X_minU(ii,4:6))*1000;
end
% Relative Distance
figure;
subplot 211; grid on; hold on;
plot(t_minU/60,position_err,'LineWidth',1)
ylabel('Relative Distance (m)')
title('Relative Distance of Chaser w.r.t Target')
% Relative Velocity
subplot 212; grid on; hold on;
plot(t_minU/60,velocity_err,'LineWidth',1)
ylabel('Relative Velocity (m/s)')
xlabel('Time (min)')
title('Relative Velocity of Chaser w.r.t Target')

% Get thrust profile and switch function
lambda_v = X_minU(:,11:13)'; 
normp = vecnorm(lambda_v);
S = normp - 1;
delta = 0.5* (1+tanh(S/rho));
figure;
subplot 211; grid on; hold on;
yyaxis left
plot(t_minU/60,delta,'DisplayName','Throttle', 'LineWidth', 1);
ylabel('Throttle (\delta)');
yyaxis right
plot(t_minU/60,S,'DisplayName','Switch Function', 'LineWidth', 1);
ylabel('Switch Function')
title('Thrust Profle and Switch Function');
legend('show','Location','best'); 
% mass
subplot 212; grid on; hold on;
plot(t_minU/60,(m0-X_minU(:,7))*1000,'LineWidth',1)
xlabel('Time (min)'); ylabel('mass (g)')
title('Mass Consumed');

% Trajectory (m)
figure; grid on; hold on; axis("equal");
plot3(X_minU(:,1)*1000,X_minU(:,2)*1000,X_minU(:,3)*1000,'-b','LineWidth',1) 
xlabel('x-position (m)'); ylabel('y-position (m)'); zlabel('z-position (m)');

end