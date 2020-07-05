clc,clear

%% Environmental Data

%Standard temperature ranges for interpolation (C)
T_wat = linspace(0,100);
T_ice = linspace(-100,0);

%Properties of water (f(T)) (SI)
rho_wat = interp1([0,20,40,60,80,100],... %Density of water
                  [999.84,998.21,992.22,983.2,971.82,958.4],...
                  T_wat); %https://www.engineersedge.com/physics/water__density_viscosity_specific_weight_13146.htm
nu_wat = 1e-6*interp1([0,20,40,60,80,100],... %Kinematic viscosity
                      [1.787,1.004,0.658,0.475,0.365,0.294],...
                      T_wat); %https://www.engineersedge.com/physics/water__density_viscosity_specific_weight_13146.htm
mu_wat = 1e-3*interp1([0,20,40,60,80,100],... %Dynamic viscosity
                      [1.793,1.002,0.6532,0.4665,0.3544,0.2818],...
                      T_wat); %https://www.engineersedge.com/physics/water__density_viscosity_specific_weight_13146.htm
cp_wat = 1e3*interp1([0,10,20,25,30,40,50,60,70,80,90,100],... %Specific heat
                    [4.2174,4.191,4.157,4.1379,4.1175,4.0737,4.0264,3.9767,3.9252,3.8729,3.8204,3.7682],...
                    T_wat);%https://www.engineeringtoolbox.com/specific-heat-capacity-water-d_660.html
k_wat = 1e-3*interp1([0,10,20,30,40,50,60,70,80,90,100],... %Thermal conductivity
                     [555.75,578.64,598.03,614.5,628.56,640.6,650.91,659.69,667.02,672.88,677.03],...
                     T_wat); %https://www.engineeringtoolbox.com/water-liquid-gas-thermal-conductivity-temperature-pressure-d_2012.html
alpha_wat = k_wat./(cp_wat.*rho_wat); %Thermal diffusivity

%Properties of ice (f(T)) (SI)
rho_ice = interp1([-100,-90,-80,-70,-60,-50,-40,-35,-30,-25,-20,-15,-10,-5,0],...
                    [925.7,924.9,924.1,923.3,922.4,921.6,920.8,920.4,920,919.6,919.4,919.4,918.9,917.5,916.2],...
                    T_ice); %https://www.engineeringtoolbox.com/ice-thermal-properties-d_576.html    
cp_ice = 1e3*interp1([-100,-90,-80,-70,-60,-50,-40,-35,-30,-25,-20,-15,-10,-5,0],...
                    [1.389,1.463,1.536,1.609,1.681,1.751,1.818,1.851,1.882,1.913,1.943,1.972,2,2.027,2.05],...
                    T_ice); %https://www.engineeringtoolbox.com/ice-thermal-properties-d_576.html
k_ice = interp1([-100,-90,-80,-70,-60,-50,-40,-35,-30,-25,-20,-15,-10,-5,0],...
                    [3.48,3.34,3.19,3.05,2.9,2.76,2.63,2.57,2.5,2.45,2.39,2.34,2.3,2.25,2.22],...
                    T_ice); %https://www.engineeringtoolbox.com/ice-thermal-properties-d_576.html
alpha_ice = k_ice./(cp_ice.*rho_ice); %Thermal diffusivity

%Constants
g_e = 1.315;  %Gravitational acceleration on Europa
g_E = 9.805;  %Gravitational acceleration on Earth
L_melt = 333.7e3; %Latent heat of fusion, a.k.a. melting enthalpy
T_melt = 0; %TODO: this should be dependent on pressure

%% Melt probe definitions
VERNE = Probe('VERNE','circle',0.35,3.71,1,448,0.12); %(Georgia Tech, SESAME, 2019 - 2021)
IceMole2 = Probe('IceMole2','square',0.15,0.7,0.3,30,0.01); %(Various, EnEx, 2014 - 2016) %TODO: Need to check this data
UlamecSample1 = Probe('UlamecSample1','circle',0.1,1,0.5,1,0.01); %(Probe for demonstrating power calculations (CG, m, dR not given), Ulamec paper on thermal probe design, 2006)
RECAS = Probe('RECAS','circle',0.16,1,0.5,1,0.01); %(Melt head tested, JLU, 2019-2020)
KowalskiSample1 = Probe('KowalskiSample1','circle',0.12,1,0.5,25,0.01); %(Probe from Kowalski 2018 paper)
%TODO: add Tunnelbot, SLUSH, Cryobot, etc.

%Test (Possibly a prototype I could build?)
m_Test = 1380*0.18*pi*((2*0.0254)^2 - (1.875*0.0254)^2) + ...
         7850*0.18*pi*(0.125*0.0254)^2 + ...
         0.11 + ...
         0.5;
Test = Probe('Test','circle',0.05,0.2,0.05,m_Test,0.05);

%% 
save("meltProbeData.mat") %Save to file