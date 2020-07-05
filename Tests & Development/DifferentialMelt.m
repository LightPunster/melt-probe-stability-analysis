addpath('..')
load('meltProbeData.mat','T_wat','T_ice','rho_wat','rho_ice','cp_wat','cp_ice','T_melt','L_melt');

shape = 'cone';
P0 = 1000;
h = 0.1;

%Functions w/ integral = 1 over cone surface area
%   3*x^2/(2*pi*h^3)

P = @(x,alpha) P0*(3*(x.^2)/(2*pi*h^3));

x = linspace(0,-h);
alpha = linspace(-pi,pi);
[X,Alpha] = meshgrid(x,alpha);
Power = P(X,Alpha);
mesh(X,Alpha,Power)


T_f = 10;
T_amb = -10;

rho_i = interp1(T_ice,rho_ice,(T_melt+T_amb)/2);
cp_i = interp1(T_ice,cp_ice,(T_melt+T_amb)/2);
rho_w = interp1(T_wat,rho_wat,(T_f+T_melt)/2);
cp_w = interp1(T_wat,cp_wat,(T_f+T_melt)/2);

P = 2000;
E = rho_i*(L_melt + cp_i*(T_melt - T_amb)) + rho_w*cp_w*(T_f - T_melt);

t = linspace(0,100,1000);
z = zeros(size(t));

for i=2:length(t)
    dt = t(i)-t(i-1);
    dDelta = (P/(2*VERNE.A_x*E))*dt;
    
    z(i) = z(i-1) + dDelta;
end

close all
figure(1)
plot(t,-z)
grid on
