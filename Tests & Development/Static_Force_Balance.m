clear, clc, close all
addpath('..')

meltProbeDataWriter
probe = VERNE;

%TODO: Include buoyancy

mu = 0.0001; %Coefficient of friction with ice wall
phi = 0*pi/180; %tile angle of melt channel w.r.t. gravity vector (rad)
theta = 2*pi/180; %tilt angle of probe w.r.t. melt channel (rad)
g = 1.315; %Gravitational acceleration (m/s^2)

theta_c = probe.ThetaContact(); %Theta to contact side walls
theta = lim(theta,[-theta_c,theta_c]); %theta may not be greater than theta_c
contact = (theta==theta_c)||(theta==-theta_c); %If theta==theta_c, probe is contacting side walls

%Notes
%F_p = reaction force at drill point
%F_w = reaction force at wall support
%par & perp refer to parallel & perpendicular (to probe axis) components of force

%Gravity
F_grav_par = -probe.m*g*cos(phi + theta);
F_grav_perp = -probe.m*g*sin(phi + theta);
F_grav = [F_grav_par; F_grav_perp];

%Moment Balance about point
%   -F_w_par*sin(gamma)*diag - F_w_perp*cos(gamma)*diag = -F_grav_perp*CG
%Force Balance
%   Parallel: F_p_par + F_w_par = -F_grav_par
%   Perpendicular: F_p_perp + F_w_perp = -F_grav_perp
%Wall Friction (assume friction at base is infinite)
%   Transform into channel(side wall) frame to relate normal & frictional forces
%   (F_w_par*cos(theta) + F_w_perp*sin(theta)) = mu*(-F_w_par*sin(theta) + F_w_perp*cos(theta))

A = [0 0        -sin(probe.gamma)*probe.diag*contact -cos(probe.gamma)*probe.diag*contact;
     1 0         1                                    0;
     0 1         0                                    1;
     0 0         cos(theta)+mu*sin(theta)             sin(theta)-mu*cos(theta)];
b = [F_grav_perp*probe.CG; -F_grav_par; -F_grav_perp; 0];
[x,status,~,~] = linEqnSol(A,b);
fprintf("%s\n",status)

if strcmp(status,'Unique solution')
    %Outputs
    F_p = [x(1); x(2)]
    F_w = [x(3); x(4)]
    F_grav

    %Balance Checks
    Moment_Balance = round(-sin(probe.gamma)*probe.diag*F_w(1) - cos(probe.gamma)*probe.diag*F_w(2) - F_grav_perp*probe.CG,5)
    Parallel_Balance = round(F_grav(1) + F_p(1) + F_w(1),5)
    Perpendicular_Balance = round(F_grav(2) + F_p(2) + F_w(2),5)
    Friction_Relation = round(F_w(1)*cos(theta) + F_w(2)*sin(theta) - mu*(-F_w(1)*sin(theta) + F_w(2)*cos(theta)),5)
end


