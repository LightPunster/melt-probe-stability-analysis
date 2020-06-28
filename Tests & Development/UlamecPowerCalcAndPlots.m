clear, clc, close all
addpath('..')

meltProbeDataWriter
%probe = UlamecSample1; V_mph = logspace(-3,1); %For verification against results in paper
probe = VERNE; V_mph = logspace(-2.3,0.9);

V = V_mph/3600;

T_amb = [50 100 150 200 250 270] - 273.15;
Lm = 334e3; %Latent heat of fusion / melting enthalpy (J/kg)
T_melt = 0; %TODO: This should change later with (pressure?)

figure
for i=1:length(T_amb)
    T = T_amb(i);
    rho = interp1(T_ice,rho_ice,(T + T_melt)/2); %From Ulamec: use averages over temp range from T_ice to T_melt
    cp = interp1(T_ice,cp_ice,(T + T_melt)/2);
    k = interp1(T_ice,k_ice,(T + T_melt)/2);

    %Ulamec 2006, pure melt probe
    a = 932; %Ws/K/m^3, fit constant
    b = 0.726;
    switch probe.cross
        case 'circle'
            A = pi*probe.r^2;
        case 'square'
            A = probe.s^2;
    end
    
    P_0 = V*(A*rho*(cp*(T_melt - T) + Lm));

    x = probe.L./(V*probe.r^2);
    if (min(x) < 5e4)
        error('Max velocity is too large for model to be valid');
    end
    if (max(x) > 1e8)
        error('Min velocity is too small for model to be valid');
    end
    P_cond = (probe.r^2)*V.*(T_melt - T)*a.*x.^b;
    P = P_0 + P_cond;
    Nu = P_0./P;
    
    subplot(2,1,1),hold on
    plot(V_mph,P/1000)
    
    subplot(2,1,2),hold on
    plot(V_mph,Nu)
end

subplot(2,1,1)
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('v(m/h)')
ylabel('P (kW)')
grid on
legend({'50 K','100 K','150 K','200 K','250 K','270 K'},'location','northwest')
title(probe.name)

subplot(2,1,2)
set(gca,'xscale','log')
xlabel('v(m/h)')
ylabel('Efficiency')
grid on
legend({'50 K','100 K','150 K','200 K','250 K','270 K'},'location','northwest')
title(probe.name)
