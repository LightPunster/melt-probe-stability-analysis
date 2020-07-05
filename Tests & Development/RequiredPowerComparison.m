clear, clc, close all
addpath('..')

meltProbeDataWriter
%probe = UlamecSample1; %For verification against results in paper
%probe = KowalskiSample1; %For verification against results in paper
probe = VERNE;

V_mph = linspace(0.1,0.7,5);
T = linspace(-20,-1,5);

[V,T] = meshgrid(V_mph, T);

[P_u,D_u,L_u,L_crit_u] = probe.RequiredPower(V,T,'Ulamec');
[P_t,D_t,L_t,L_crit_t] = probe.RequiredPower(V,T,'Talalay');
[P_k,D_k,L_k,L_crit_k] = probe.RequiredPower(V,T,'Kowalski');
[P_l,D_l,L_l,L_crit_l] = probe.RequiredPower(V,T,'Li');

%%
figure(1),hold on
mesh(V,T,P_u,'EdgeColor','r')
mesh(V,T,P_t,'EdgeColor','g')
mesh(V,T,P_k,'EdgeColor','b')
mesh(V,T,P_l,'EdgeColor','k')
%set(gca,'xscale','log')%, set(gca,'zscale','log')
grid on
xlabel('V (m/hr)'), ylabel('T (C)'), zlabel('P (W)')
legend('Ulamec','Talalay','Kowalski','Li')
title('Required Power')

figure(2),hold on
mesh(V,T,D_u,'EdgeColor','r')
mesh(V,T,D_t,'EdgeColor','g')
mesh(V,T,D_k,'EdgeColor','b')
mesh(V,T,D_l,'EdgeColor','k')
%set(gca,'xscale','log')%, set(gca,'zscale','log')
grid on
xlabel('V (m/hr)'), ylabel('T (C)'), zlabel('P (W)')
legend('Ulamec','Talalay','Kowalski','Li')
title('Descent Power')

figure(3),hold on
mesh(V,T,L_u,'EdgeColor','r')
mesh(V,T,L_t,'EdgeColor','g')
mesh(V,T,L_k,'EdgeColor','b')
mesh(V,T,L_l,'EdgeColor','k')
%set(gca,'xscale','log')%, set(gca,'zscale','log')
grid on
xlabel('V (m/hr)'), ylabel('T (C)'), zlabel('P (W)')
legend('Ulamec','Talalay','Kowalski','Li')
title('Lateral Power')

figure(4),hold on
mesh(V,T,L_crit_u,'EdgeColor','r')
mesh(V,T,L_crit_t,'EdgeColor','g')
mesh(V,T,L_crit_k,'EdgeColor','b')
mesh(V,T,L_crit_l,'EdgeColor','k')
%set(gca,'xscale','log')%, set(gca,'zscale','log')
grid on
xlabel('V (m/hr)'), ylabel('T (C)'), zlabel('L_{crit} (m)')
legend('Ulamec','Talalay','Kowalski','Li')
title('Critical Refreezing Length')
xlimits = xlim;
ylimits = ylim;
axis([xlimits ylimits 0 5*probe.L])

