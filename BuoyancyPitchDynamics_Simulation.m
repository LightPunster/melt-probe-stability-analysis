%% Setup
meltProbeDataWriter
clc, clear, close all

probe = load('meltProbeData.mat',"VERNE").VERNE;

%linearize = true;
linearize = false;

%Temperature (C)
    T = 5; 

%Time Range (s)
    t = linspace(0,100,10000);

%Initial state
theta0 = 1*pi/180;
omega0 = 0;
x0 = [theta0,omega0];
    
%Disturbance
    %Force at CG (N)
        %F_dist = 0;
        F_dist = 0;
        %F_dist = 0.6*(rand(1,length(t)) - 0.5); %White noise, will be important for Monte Carlos later
        %F_dist = 3*sin(2*t); %(this is highly attenuated)
        %F_dist = 3*sin(t/2); %(this is close to resonance)
    %Time of disturbance (s)
        t_dist = 40;
        Tau = F_dist.*(t<=t_dist).*probe.CG;

%% Performance Prediction
n = 100;

CGs = linspace(probe.L/10,9*probe.L/10,n);
ms = linspace(0.5*probe.m,1.5*probe.m,n);
Ls = linspace(0.5*probe.L,1.5*probe.L,n);
switch probe.cross
    case('circle'), ws = linspace(0.8*probe.d,1.2*probe.d,n);
    case('square'), ws = linspace(0.8*probe.s,1.2*probe.s,n);
    otherwise, error('Invalid cross section specified.');
end
EigenVals_CG = zeros(2,length(CGs));
EigenVals_m = zeros(2,length(ms));
EigenVals_L = zeros(2,length(Ls));
EigenVals_w = zeros(2,length(ws));

probe0 = probe; %Save original probe
for i=1:n
    probe.CG = CGs(i); probe = probe.Update(); %Update parameters
    [~,~,EigenVals_CG(:,i)] = probe.WaterPocketPitchDynamics_Linear(T);
end
probe = probe0; %Reload original probe
for i=1:n
    probe.m = ms(i); probe = probe.Update(); %Update parameters
    [~,~,EigenVals_m(:,i)] = probe.WaterPocketPitchDynamics_Linear(T);
end
probe = probe0; %Reload original probe
for i=1:n
    probe.L = Ls(i); probe = probe.Update(); %Update parameters
    [~,~,EigenVals_L(:,i)] = probe.WaterPocketPitchDynamics_Linear(T);
end
probe = probe0; %Reload original probe
for i=1:n
    probe.d = ws(i); probe = probe.Update(); %Update parameters
    [~,~,EigenVals_w(:,i)] = probe.WaterPocketPitchDynamics_Linear(T);
end
probe = probe0; %Reload original probe

%% Simulation
%tic %Track simulation time
if linearize
    [A,B,~] = probe.WaterPocketPitchDynamics_Linear(T);
    sys = ss(A,B,[],[]);
    [~,~,x] = lsim(sys,Tau,t,x0,'foh');
else
    [~,~,x] = nlsim(@(x,u) probe.WaterPocketPitchDynamics_Nonlinear(x,u,T),Tau,t,x0);
end
%toc %Report simulation time

%% Plot performance prediction
colors = ['r','g','b','c','m','y','k'];
figure

subplot(2,2,1),hold on
for i=1:size(EigenVals_CG,1)
    plot(CGs,real(EigenVals_CG(i,:)),colors(i))
end
vline(probe.CG,'k--'), vline(probe.CB,'b--'), hline(0,'k')
xlabel('CB (m)'), ylabel('Re(\lambda)'),grid on
legend('\lambda_{x1}','\lambda_{x2}','CG','CB','Marginal Stability')

subplot(2,2,2),hold on
for i=1:size(EigenVals_m,1)
    plot(ms,real(EigenVals_m(i,:)),colors(i))
end
vline(probe.m,'k--'), hline(0,'k')
xlabel('m (kg)'), ylabel('Re(\lambda)'),grid on
legend('\lambda_{x1}','\lambda_{x2}','m','Marginal Stability')

subplot(2,2,3),hold on
for i=1:size(EigenVals_L,1)
    plot(Ls,real(EigenVals_L(i,:)),colors(i))
end
vline(probe.L,'k--'), hline(0,'k')
xlabel('L (m)'), ylabel('Re(\lambda)'),grid on
legend('\lambda_{x1}','\lambda_{x2}','L','Marginal Stability')

subplot(2,2,4),hold on
for i=1:size(EigenVals_w,1)
    plot(ws,real(EigenVals_w(i,:)),colors(i))
end
ylabel('Re(\lambda)'),grid on
switch probe.cross
    case('circle')
        vline(probe.d,'k--'), hline(0,'k')
        xlabel('d (m)')
        legend('\lambda_{x1}','\lambda_{x2}','d','Marginal Stability')
    case('square')
        vline(probe.s,'k--'), hline(0,'k')
        xlabel('s (m)')
        legend('\lambda_{x1}','\lambda_{x2}','s','Marginal Stability')
    otherwise, error('Invalid cross section specified.');
end

%% Plot simulation
theta = x(:,1); %Angle from verticle
theta_d = x(:,2); %Angular velocity
theta_c = probe.ThetaContact(); %Pitch at which side walls are contacted

%End plotted output when probe hits side walls of melt channel
i_UL = find(theta>min(theta_c),1,'first');
i_LL = find(theta<-min(theta_c),1,'first');
i_end = min([i_UL,i_LL,length(t)]);

theta = theta(1:i_end);
theta_d = theta_d(1:i_end);
t = t(1:i_end);
Tau = Tau(1:i_end);

figure
subplot(3,1,1)
plot(t,Tau),xlabel('t (s)'),ylabel('T_{dist} (N-m)'),grid on
subplot(3,1,2)
plot(t,(180/pi)*theta),xlabel('t (s)'),ylabel('\theta (deg)'),grid on
subplot(3,1,3)
plot(t,(180/pi)*theta_d),xlabel('t (s)'),ylabel('\theta_d (deg/s)'),grid on

%% Calculate Reynolds (verifies pressure drag is insignificant)
figure
[R,Theta_d,L] = RotatingReynoldsNumber(probe,theta_d,T);
mesh(L,Theta_d,R)

%%
DragForceEstimates