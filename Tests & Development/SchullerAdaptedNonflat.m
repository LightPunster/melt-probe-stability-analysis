addpath('..')
meltProbeDataWriter
probe = VERNE;

T_amb = -10;
T_nose = 20;
V = (1/3600)*[0.57,0.2,0];
d0 = delta(T_nose,T_amb,norm(V)); %Sanity check

switch (probe.cross)
    case 'circle'
        dim = probe.r;
    case 'square'
        dim = probe.s;
end
y = linspace(-dim,dim,20);
z = y;
[Y,Z] = meshgrid(y,z);
X_p = plate(Y);
X_d = dome(dim,Y,Z);
X_c = cone(dim,dim,Y,Z);
switch(probe.cross)
    case 'circle'
        X_p((Y.^2 + Z.^2)>probe.r^2)=NaN;
        X_d((Y.^2 + Z.^2)>probe.r^2)=NaN;
        X_c((Y.^2 + Z.^2)>probe.r^2)=NaN;
end

close all

figure(1)
mesh(X_p,Y,Z)
for i=1:numel(X_p)
    v_x = proj(V,plateNorm());
    plotVect(100*(v_x/norm(v_x))*delta(T_nose,T_amb,norm(v_x)),...
        [X_p(i),Y(i),Z(i)],'k')
end
xlabel('x'),ylabel('y'),zlabel('z')
axis([-2*dim,2*dim,-2*dim,2*dim,-2*dim,2*dim])

figure(2)
mesh(X_d,Y,Z)
for i=1:numel(X_d)
    v_x = proj(V,domeNorm(probe.r,X_d(i),Y(i),Z(i)));
    plotVect(100*(v_x/norm(v_x))*delta(T_nose,T_amb,norm(v_x)),...
        [X_d(i),Y(i),Z(i)],'k')
end
xlabel('x'),ylabel('y'),zlabel('z')
axis([-2*dim,2*dim,-2*dim,2*dim,-2*dim,2*dim])

figure(3)
mesh(X_c,Y,Z)
for i=1:numel(X_c)
    v_x = proj(V,coneNorm(probe.r,probe.r,Y(i),Z(i)));
    plotVect(100*(v_x/norm(v_x))*delta(T_nose,T_amb,norm(v_x)),...
        [X_c(i),Y(i),Z(i)],'k')
end
xlabel('x'),ylabel('y'),zlabel('z')
axis([-2*dim,2*dim,-2*dim,2*dim,-2*dim,2*dim])


function X = plate(Y)
    X = zeros(size(Y));
end

function n = plateNorm()
    n = [1,0,0];
end

function X = dome(r,Y,Z)
    X = real(sqrt(r.^2 - Y.^2 - Z.^2)) - r;
end

function n = domeNorm(r,x,y,z)
    n = [(x+r)/r,y/r,z/r];
end

function X = cone(r,l,Y,Z)
    X = -(l/r)*sqrt(Y.^2 + Z.^2);
end

function n = coneNorm(r,l,y,z)
    n = (1/sqrt(1+l/r))*[1, sqrt(l/r)*sin(atan2(y,z)), sqrt(l/r)*cos(atan2(y,z))];
end

function p = proj(a,b) %Projection of a onto b
    p = (dot(a,b)/norm(b)^2)*b;
end

%From Schuller 2016 - Curvilinear melting ...
function d = delta(T_n,T_a,V)
    load('meltProbeData.mat','T_wat','T_ice','T_melt','k_wat','cp_ice','rho_ice','L_melt')
    k_w = interp1(T_wat,k_wat,(T_melt+T_n)/2);
    cp_i = interp1(T_ice,cp_ice,(T_melt+T_a)/2);
    rho_i = interp1(T_ice,rho_ice,(T_melt+T_a)/2);
    
    d = k_w*(T_n - T_melt)/(rho_i*(L_melt + cp_i*(T_melt - T_a))*V);
end
    
