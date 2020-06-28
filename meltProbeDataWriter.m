clc,clear

%Environmental Data
T_H2O = linspace(0,100); %Standard temperature range for interpolation
T_ice = linspace(-100,0);
%https://www.engineersedge.com/physics/water__density_viscosity_specific_weight_13146.htm
rho_H2O = interp1([0,20,40,60,80,100],[999.84,998.21,992.22,983.2,971.82,958.4],T_H2O); %density of water
nu_H2O = 1e-6*interp1([0,20,40,60,80,100],[1.787,1.004,0.658,0.475,0.365,0.294],T_H2O); %kinematic viscosity of water
mu_H2O = 1e-3*interp1([0,20,40,60,80,100],[1.793,1.002,0.6532,0.4665,0.3544,0.2818],T_H2O); %dynamic viscosity of water
g_e = 1.315;  %gravitational acceleration on Europa
L_melt = 334e3; %Latent heat of fusion / melting enthalpy (J/kg)
T_melt = 0;

%https://www.engineeringtoolbox.com/ice-thermal-properties-d_576.html
cp_ice = 1e3*interp1([-100,-90,-80,-70,-60,-50,-40,-35,-30,-25,-20,-15,-10,-5,0],...
                     [1.389,1.463,1.536,1.609,1.681,1.751,1.818,1.851,1.882,1.913,1.943,1.972,2,2.027,2.05],...
                     T_ice);
k_ice = interp1([-100,-90,-80,-70,-60,-50,-40,-35,-30,-25,-20,-15,-10,-5,0],...
                     [3.48,3.34,3.19,3.05,2.9,2.76,2.63,2.57,2.5,2.45,2.39,2.34,2.3,2.25,2.22],...
                     T_ice);
rho_ice = interp1([-100,-90,-80,-70,-60,-50,-40,-35,-30,-25,-20,-15,-10,-5,0],...
                     [925.7,924.9,924.1,923.3,922.4,921.6,920.8,920.4,920,919.6,919.4,919.4,918.9,917.5,916.2],...
                     T_ice);    

VERNE = Probe('VERNE','circle',0.35,3.71,1,448,0.12); %(Georgia Tech, SESAME, 2019 - 2021)
IceMole2 = Probe('IceMole2','square',0.15,0.7,0.3,30,0.01); %(Various, EnEx, 2014 - 2016) %TODO: Need to check this data
UlamecSample1 = Probe('UlamecSample1','circle',0.1,1,0.5,1,0.01); %(Probe for demonstrating power calculations (CG, m, dR not given), Ulamec paper on thermal probe design, 2006)
RECAS = Probe('RECAS','circle',0.16,1,0.5,1,0.01); %(Melt head tested, JLU, 2019-2020)
%TODO: add Tunnelbot, SLUSH, Cryobot, etc.

%Test (Possibly a prototype I could build?)
m = 1380*0.18*pi*((2*0.0254)^2 - (1.875*0.0254)^2) + ...
         7850*0.18*pi*(0.125*0.0254)^2 + ...
         0.11 + ...
         0.5;
Test = Probe('Test','circle',0.05,0.2,0.05,m,0.05);

save("meltProbeData.mat") %Save to file