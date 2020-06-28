load('meltProbeData.mat','nu_H2O','rho_H2O','T_H2O','g_e');

T = 5;
nu = interp1(T_H2O,nu_H2O,T);
rho = interp1(T_H2O,rho_H2O,T);

%theta_d = logspace(-7,-1);
T_d = zeros(size(theta_d));

for i=1:length(theta_d)
    T_d(i) = probe.IntegratedPressureDragMoment(rho,nu,theta_d(i),1);
end

theta_d2 = (theta_d.^2)./sign(theta_d);
B2 = max(T_d)/max(abs(theta_d2));
T_d_fit = -B2*theta_d2;
fprintf('%.10f\n',B2)
perc_error = 100*abs((T_d - T_d_fit)./T_d);
figure
histogram(perc_error)

figure
subplot(2,1,1)
plot(abs(theta_d2),abs(T_d),abs(theta_d2),abs(T_d_fit),'--');
set(gca,'xscale','log')
set(gca,'yscale','log')
grid on

subplot(2,1,2),hold on
plot(abs(theta_d2),RotatingReynoldsNumber(probe,theta_d,T,probe.L/4))
plot(abs(theta_d2),RotatingReynoldsNumber(probe,theta_d,T,probe.L/2))
plot(abs(theta_d2),RotatingReynoldsNumber(probe,theta_d,T,3*probe.L/4))
plot(abs(theta_d2),RotatingReynoldsNumber(probe,theta_d,T,probe.L))
set(gca,'xscale','log')
set(gca,'yscale','log')
grid on