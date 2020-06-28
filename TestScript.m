clc,clear,close all

load 'meltProbeData.mat'
probe = VERNE;

theta_c = probe.ThetaContact(); %Theta to contact side walls
mu = 0.001;
rho = 997;
g = 1.315;
B = 200;
Beta = 0;
D = 10000;
T = 5;

Phi = (-3:0.5:3)*pi/180;
Theta = zeros(1,length(Phi));
Lambda_max = zeros(1,length(Phi));
F_g = zeros(2,length(Phi)); %Gravity
T_g = zeros(1,length(Phi));
F_b = zeros(2,length(Phi)); %Buoyancy
T_b = zeros(2,length(Phi));
F_p = zeros(2,length(Phi)); %Contact at point
F_w = zeros(2,length(Phi)); %Contact at wall

figure(1),hold on
for i=1:length(Phi)
    [~,~,Lambda] = probe.WaterPocketPitchDynamics_Linear(T);
    Lambda_max(i) = max(real(Lambda));
    if Lambda_max(i)<0 %Pitch stability check
        Theta(i) = -Phi(i); %If stable, probe will remain vertical
    else
        Theta(i) = theta_c; %If stable, probe will remain vertical
    end
    Theta(i) = lim(Theta(i),[-theta_c,theta_c]); %theta may not be greater than theta_c

    F_g(:,i) = probe.m*g*[0;-1];
    F_b(:,i) = BuoyantForce_2D(probe,rho,g,Beta,D,Theta(i)+Phi(i));
    [F_p(:,i),F_w(:,i),~] = StaticContactForces_2D(probe,Theta(i),Phi(i),mu,F_g(:,i)+F_b(:,i));
    
end
%[Phi*180/pi; Theta*180/pi; (Phi+Theta)*180/pi; Lambda_max; F_g; F_b; F_p; F_w]

figure(1)
subplot(2,2,1),hold on
    for i=1:length(Phi)
        plot([0,F_g(1,i)],[0,F_g(2,i)])
        try
            legappend(['\phi=' num2str(Phi(i)*180/pi)])
        catch
            legend(['\phi=' num2str(Phi(i)*180/pi)],'location','west')
        end
    end
    title('F_g (N)'),xlabel('x'),ylabel('y'),grid on
subplot(2,2,2),hold on
    for i=1:length(Phi)
        plot([0,F_b(1,i)],[0,F_b(2,i)])
    end
    title('F_b (N)'),xlabel('x'),ylabel('y'),grid on
subplot(2,2,3),hold on
    for i=1:length(Phi)
        plot([0,F_p(1,i)],[0,F_p(2,i)])
    end
    title('F_p (N)'),xlabel('x'),ylabel('y'),grid on
subplot(2,2,4),hold on
    for i=1:length(Phi)
        plot([0,F_w(1,i)],[0,F_w(2,i)])
    end
    title('F_w (N)'),xlabel('x'),ylabel('y'),grid on