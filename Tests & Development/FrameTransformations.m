clc, clear, close all

r = 0.175;
cross = 'circle';
shape = 'dome'; l = 2*r;
theta = 1;
h0 = 0.01;

n = 100;

y = linspace(-1.2*r,1.2*r,n);
z_c = linspace(-1.2*r,1.2*r,n);
[Y,Z_c] = meshgrid(y,z_c);
Z_p = Z_c;
switch shape
    case 'plate'
        B_m_c = plate(r,Y,Z_c);%Probe surface in probe frame is same as melt boundary in channel frame
        S_p_p = plate(r,Y,Z_p); %probe surface in probe frame
    case 'dome'
        B_m_c = dome(r,Y,Z_c);%Probe surface in probe frame is same as melt boundary in channel frame
        S_p_p = dome(r,Y,Z_p); %probe surface in probe frame
    case 'cone'
        B_m_c = cone(r,l,Y,Z_c);%Probe surface in probe frame is same as melt boundary in channel frame
        S_p_p = cone(r,l,Y,Z_p); %probe surface in probe frame
end
%Use transformation to calculate probe surface in channel frame
S_p_c = zeros(size(S_p_p)); 
Z_p_c = zeros(size(Z_p));
T = PtoC(theta*pi/180,-h0);
for i=1:numel(S_p_p)
    v_p = T*[S_p_p(i); Y(i); Z_p(i); 1];
    S_p_c(i) = v_p(1);
    Z_p_c(i) = v_p(3);
end

S_p_c2 = zeros(size(S_p_c));
for i=1:n
    S_p_c2(:,i) = interp1(Z_p_c(:,i),S_p_c(:,i),Z_c(:,i));
end

H = max(0,B_m_c - S_p_c2);

switch(cross)
    case 'circle'
        S_p_c2((Y.^2 + Z_c.^2)>r^2)=NaN;
        B_m_c((Y.^2 + Z_c.^2)>r^2)=NaN;
        H((Y.^2 + Z_c.^2)>r^2)=NaN;
end

figure(1),hold on
mesh(Y,Z_c,B_m_c)
mesh(Y,Z_c,S_p_c2)
axis equal
grid on

figure(2),hold on
mesh(Y,Z_c,H)
grid on
title('Melt Film Thickness')

function X = plate(r,Y,Z)
    X = zeros(size(Y));
end

function X = dome(r,Y,Z)
    X = real(sqrt(r.^2 - Y.^2 - Z.^2) - r);
end

function X = cone(r,l,Y,Z)
    X = -(l/r)*sqrt(Y.^2 + Z.^2);
end

function T_y = RotY(theta)
    T_y = [ cos(theta) 0 -sin(theta)  0
            0          1  0           0
            sin(theta) 0  cos(theta)  0
            0          0  0           1 ];
end

function T_t = TransX(x)
    T_t = [ 1 0 0 x
            0 1 0 0
            0 0 1 0
            0 0 0 1 ];
end

function T_t = UnifScale(k)
    T_t = [ k 0 0 0
            0 k 0 0
            0 0 k 0
            0 0 0 1 ];
end

function T_cp = CtoP(theta,h0)
    T_cp = TransX(-h0)*RotY(theta);
end
function T_pc = PtoC(theta,h0)
    T_pc = TransX(h0)*RotY(-theta);
end