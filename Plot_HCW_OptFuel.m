function Plot_HCW_OptFuel
%%  < File Description >
%    Author:        Ruthvik Bommena
%    File Name:     Plot_HCW_OptFuel.m
%    Compiler:      MATLAB R2022b
%    Date:          22 March, 2023
%    Affiliation:   Department of Aerospace Engineering, University of Illinois Urbana-Champaign.
%    Description:   Function plots HCW optimal fuel transfer for chaser to rendezvous with the target.

close all; clear; clc;

% Parameters
aT = 42165; % km (Target orbit radius - GEO)
mu = 398600; % km^3/s^2
m0 = 50; % kg
T =  0.0015/1000; % kN (MiXI Thruster 1.5mN)
Isp = 4190; % s
g0 = 9.8/1000; % km/s^2
c = Isp*g0; % km/s

% Boundary conditions (km and km/s)
x0 = [0 1 0 0 0 0]';
xf = [0 0 0 0 0 0]';

% Transfer time
tf = 200; % min
tf = tf*60; % sec

% OptFuel transfer
[t_minU,X_minU,lam0,rho] = Solve_HCWOptFuel(x0,xf,tf,m0,mu,T,c,aT);

% Plots
plots(t_minU,X_minU,rho,m0,c,x0,xf)

end


%% Function Plots
function plots(t_minU,X_minU,rho,m0,c,x0,xf)

x0(1:3) = x0(1:3)*1000; % change to meters
xf(1:3) = xf(1:3)*1000; % change to meters

% Get relative distance and velocity of chaser w.r.t target
% Get switch function and thrust profile
for ii = 1:length(X_minU)
    rel_position(ii) = norm(X_minU(ii,1:3))*1000; % m
    rel_velocity(ii) = norm(X_minU(ii,4:6))*1000; % m/s

    S(ii) = norm(X_minU(ii,11:13))*c/X_minU(ii,7) + X_minU(ii,14) - 1; % switch function
    delta(ii) = 0.5* (1+tanh(S(ii)/rho)); % throttle
end

% Relative Distance
figure;
subplot 211; grid on; hold on;
plot(t_minU/60,rel_position,'LineWidth',1)
ylabel('Relative Distance (m)'); xlabel('Time (min)');
title('Relative Distance of Chaser w.r.t Target')
% Relative Velocity
subplot 212; grid on; hold on;
plot(t_minU/60,rel_velocity,'LineWidth',1)
ylabel('Relative Velocity (m/s)')
xlabel('Time (min)')
title('Relative Velocity of Chaser w.r.t Target')

% thrust profile and switch function
figure;
subplot 211; grid on; hold on;
yyaxis left
plot(t_minU/60,delta,'DisplayName','Throttle', 'LineWidth', 1);
ylabel('Throttle (\delta)');
yyaxis right
plot(t_minU/60,S,'DisplayName','Switch Function', 'LineWidth', 1);
ylabel('Switch Function'); xlabel('Time (min)');
title('Thrust Profle and Switch Function');
legend('show','Location','best'); 
% mass
subplot 212; grid on; hold on;
plot(t_minU/60,(m0-X_minU(:,7))*1000,'LineWidth',1)
xlabel('Time (min)'); ylabel('mass (g)')
title('Mass Consumed');

% Trajectory (m)
figure; grid on; hold on; axis("equal");
plot3(x0(1), x0(2), x0(3), '*k', 'LineWidth', 1.5, 'Markersize', 15, 'DisplayName', 'Initial Position');
plot3(xf(1), xf(2), xf(3), '*r', 'LineWidth', 1.5, 'Markersize', 15, 'DisplayName', 'Final Position');
plot3(X_minU(:,1)*1000,X_minU(:,2)*1000,X_minU(:,3)*1000,'-b','LineWidth',2,'DisplayName', 'Trajectory') 
xlabel('x-position (m)'); ylabel('y-position (m)'); zlabel('z-position (m)');
legend('show','Location','best'); 

end