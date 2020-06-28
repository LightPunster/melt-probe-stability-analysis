clear,clc,close all

param.r = 0.175;
param.L = linspace(2,3,100);
param.R_delta = 0.12;
param.theta_c = 0; % As the output, inital value doesn't matter.

param = multivariateTrade(param, @Residual, 'theta_c');

% gamma = (180/pi)*atan(param.r./param.L);
% figure,plot(param.r,gamma)
% 
% R = param.r + param.R_delta;
% psi = (180/pi)*acos(R./sqrt(param.r.^2 + param.L.^2));
% figure,plot(param.r,psi)

%%
function res = Residual(x,p,o)
  p.(o) = x;
  
  R = p.r + p.R_delta;
  gamma = atan(p.r./p.L);
  psi = acos(R./sqrt(p.r.^2 + p.L.^2));
  res = pi/2 - p.theta_c*(pi/180) - gamma - psi;
end
