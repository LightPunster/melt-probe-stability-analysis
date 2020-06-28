function [R,Theta_d,L] = RotatingReynoldsNumber(probe,theta_d,T,varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    
    
    load('meltProbeData.mat','nu_H2O','T_H2O');
    nu = interp1(T_H2O,nu_H2O,T);
    
    %Section on flow around airfoils in: https://en.wikipedia.org/wiki/Reynolds_number
    c = probe.d;
    switch(nargin)
        case 3
            L = linspace(0,probe.L);
            [L,Theta_d] = meshgrid(L,theta_d);
            V = abs(Theta_d).*L;
            R = V*c/nu;
        case 4
%             L = linspace(0,probe.L);
%             [L,Omega] = meshgrid(L,theta_d);
%             V = abs(Omega).*L;
%             R = V*c/nu;
%             if strcmp(varargin{1},'plot')
%                 mesh(L,Omega,R)
%                 xlabel('l (m)'),ylabel('\theta_d (rad/s)'),zlabel('R')
%                 %title('Reynolds Number')
%                 %set(gca,'zscale','log')
%             end
%         case 5
%             V = abs(theta_d)*varargin{2};
            Theta_d = theta_d;
            L = varargin{1};
            V = abs(Theta_d).*L;
            R = V*c/nu;
%             if strcmp(varargin{1},'plot')
%                 plot(theta_d,Re)
%                 xlabel('\theta_d (rad/s)'),ylabel('Re')
%                 %title('Reynolds Number')
%                 %set(gca,'zscale','log')
%             end 
        otherwise
            error('Expected 2 or 3 inputs, got %d.',nargin)
            %error('Expected 2, 3, or 4 inputs, got %d.',nargin)
    end
end

