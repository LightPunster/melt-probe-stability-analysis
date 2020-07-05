
%TODO: Account for Beta factor at melt head contact

classdef Probe
    %Probe stores the physical properties of a long, thin ice melt probe
    %which are needed to perform dynamic stability analysis
    %
    %   Assumptions:
    %       1. Constant cross section from nose to tail
    %       2. Center of buoyancy is at center of volume 
    %       3. Width is small relative to length
    %       4. Probe is symmetric about longitudinal axis
    %       5. Melt channel has same cross-section shape as melt probe
    %
    %   Assessment of assumptions
    %           The deviation is most likely to be seen at the nose (which,
    %       for example, may be a cone, parabaloid, ellipsoid, or similar).
    %       Since the nose will makeup a very small portion of the probe's
    %       total surface area and an even smaller portion of its volume,
    %       the impact of assuming a constant cross section is minimal. It
    %       will manifest in the model as a slightly lower-than-true center
    %       of buoyancy (lower stiffness), a slightly lower-than-true
    %       center of hydrodynamic forces like drag (lower damping) (this
    %       should be almost insignificant compared to other forces,
    %       regardless), a slightly higher than true moment of inertia, and
    %       a slightly higher than true Added Mass tensor. The first two
    %       effects will lead to a conservative estimate of stability which
    %       may be adjusted for, if desired, by altering the initialized
    %       length and CG of the probe. The second two effects should be
    %       very small because of the melt head's proximity to the axis of
    %       rotation. An exception to this may occur if the melt head is an
    %       extremely large portion of the probe's mass, but this seems
    %       unlikely from examining currently existing designs,
    %       particularly given the assumption #3. Regardless, a
    %       deviation in moment of inertia will not affect the stability,
    %       but may have a small impact on the magnitude of angular
    %       acceleration.
    %
    %   Author: Nathan Daniel
    %   Date: 06/26/2020
    
    properties
        %Initialized
        name %Name of melt probe
        cross %Constant cross-sectional shape
        d = 'N/A'; %Diameter of circular cross-section probe (m)
        s = 'N/A'; %Side width for square cross-section probe (m)
        L %Total length (m)
        CG %Center of gravity (m from nose)
        m %Total mass of probe (kg)
        R_delta %TODO: Parameter which may be later replaced, but described distance from side of probe to walls of melt channel, in m
        
        %Calculated by Update
        r = 'N/A'; %Radius of circular cross-section probe (m)
        a = 'N/A'; %1/2 side width for square cross-section probe (m)
        CB %Center of buoyancy (m from nose)
        AR %Aspect ratio
        A_x %Horizontal cross-sectional area (m^2)
        A_y %Vertical cross-sectional area (m^2)
        V %Volume (m^3)
        diag %Length of diagonal (line from drill point to upper edge) (m)
        gamma %Angle from longitudinal axis to diagonal (radians) %TODO: Maybe throw warning if this is too large?
        I %3x3 Inertia tensor (kg*m^2)
        B1_mu %Damping coefficient(s), normalized by TODO (dynamic, kinematic?) viscosity of fluid TODO (units?)
        B2_rho %Drag coefficient(s), normalized by density of fluid (N/(m/s)^2/(kg/m^3))
        M_rho %6x6 Added mass tensor, normalized by density of fluid (kg, m^2, or kg*m^2 depending on entry)
    end
    
    methods
        function obj = Probe(name,cross,w,L,CG,m,R_delta)
            %Probe Construct an instance of the probe class
            %   Parameters: 
            %       name = name of melt probe
            %       cross = constant cross-sectional shape
            %       w = diameter (m) (for cross='circle') or side width (m)
            %           (for cross = 'square')
            %       L = total length (m)
            %       CG = center of gravity (m from nose)
            %   Author: Nathan Daniel
            %   Date: 06/26/2020
    
            obj.name = name;
            obj.cross = cross;
            switch cross
                case 'circle', obj.d = w;
                case 'square', obj.s = w;
                otherwise, error("%s is not a valid value for 2nd constructor parameter 'cross'\n...(must be 'circle' or 'square')",cross);
            end
            obj.L = L;
            obj.CG = CG;
            obj.m = m;
            obj.R_delta = R_delta;
            obj = obj.Update();
        end
        
        function obj = Update(obj)
            %Update calculates derived parameters of the Probe based on
            %the lower-level ones defined in the Probe constructor.
            %   Detailed explanation TODO
            %
            %   Author: Nathan Daniel
            %   Date: 06/26/2020
            
            switch(obj.cross)
                case 'circle'
                    if ischar(obj.d) %Switching from square to circle
                        obj.d = obj.s; obj.s = 'N/A';
                    end
                        
                    %Dimensional
                    obj.r = obj.d/2;
                    obj.CB = obj.L/2;
                    obj.AR = obj.d/obj.L;
                    obj.A_x = pi*(obj.r^2);
                    obj.A_y = obj.d*obj.L;
                    obj.V = obj.A_x*obj.L;
                    
                    %For calculating contact pitch
                    obj.diag = sqrt(obj.r^2 + obj.L^2);
                    obj.gamma = acos(obj.L/obj.diag);

                    %Inertia Tensor Reference: https://en.wikipedia.org/wiki/List_of_moments_of_inertia#List_of_3D_inertia_tensors
                    I_long = 0.5*obj.m*obj.r^2; %Moment of inertia about longitudinal axis
                    I_lat = (1/12)*(0.5*obj.m)*(3*obj.r^2 + obj.CG^2) + ... %Moment of inertia of mass below CG about drill point
                            (1/12)*(0.5*obj.m)*(3*obj.r^2 + (obj.L-obj.CG)^2) + (0.5*obj.m)*obj.CG^2; %Moment of inertia of mass above CG about drill point (via parallel axis theorem)
                    obj.I = [I_long 0     0
                             0      I_lat 0
                             0      0     I_lat ];

                case 'square'
                    if ischar(obj.s) %Switching from circle to square
                        obj.s = obj.d; obj.d = 'N/A';
                    end
                    
                    %Dimensional
                    obj.a = obj.s/2;
                            obj.CB = obj.L/2;
                    obj.AR = obj.s/obj.L;
                    obj.A_x = obj.s^2;
                    obj.A_y = [ obj.s*obj.L %Min cross section (rectangular)
                                sqrt(2)*obj.s*obj.L ]; %Max (diamond)
                    obj.V = obj.A_x*obj.L;

                    %For calculating contact pitch (TODO)
                    obj.diag = [sqrt(obj.a^2 + obj.L^2) %First entry is edge diagonal
                                sqrt(2*obj.a^2 + obj.L^2) ]; %Second entry is corner diagonal
                    obj.gamma = acos(obj.L./obj.diag); %First entry is edge, second is corner
                    
                    %Inertia Tensor Reference: https://en.wikipedia.org/wiki/List_of_moments_of_inertia#List_of_3D_inertia_tensors
                    I_long = (1/12)*obj.m*(obj.s^2 + obj.s^2); %Moment of inertia about longitudinal axis
                    I_lat = (1/12)*(0.5*obj.m)*(obj.s^2 + obj.CG^2) + (0.5*obj.m)*(0.5*obj.CG)^2 + ... %Moment of inertia of mass below CG about drill point (via parallel axis theorem)
                            (1/12)*(0.5*obj.m)*(obj.s^2 + (obj.L-obj.CG)^2) + (0.5*obj.m)*(obj.CG + (obj.L-obj.CG)/2)^2; %Moment of inertia of mass above CG about drill point (via parallel axis theorem)
                    obj.I = [I_long 0     0
                             0      I_lat 0
                             0      0     I_lat ];
                otherwise
                    error("%s is not a valid value for property 'cross'\n...(must be 'circle' or 'square')",obj.cross);
            end
            %Hydrodynamic Force Coefficients
            [obj.B1_mu,obj.B2_rho] = RotationalDragCoeffs(obj);
            obj.M_rho = obj.AddedMassTensor();
        end
        
        function M_over_rho = AddedMassTensor(obj) %TODO: Add help
            %UNTITLED3 Summary of this function goes here
            %   Detailed explanation goes here

            %Reference: https://ocw.mit.edu/courses/mechanical-engineering/2-016-hydrodynamics-13-012-fall-2005/readings/2005reading6.pdf

            %Coefficients for 2D cross sections: Page 4
            %NOTE: normalized by (constant) density rho
            switch(obj.cross)
                case 'circle'
                    A = [obj.d^2 0         0
                         0         obj.d^2 0
                         0         0         0 ];
                case 'square'
                    A = [1.51*pi*obj.a^2 0                 0
                         0                 1.51*pi*obj.a^2 0
                         0                 0                 0.234*pi*obj.a^4 ];
            end

            %2D Cross-section --> 3D: Page 9
            %NOTE: interpretation below valid only for constant cross
            %section (Assumption 1)
            %NOTE: may also only be valid for long/narrow, because it
            %neglects effects of top and bottom surfaces
            m_22 = A(1,1)*obj.L;
            m_23 = -A(1,2)*obj.L;
            m_24 = A(1,3)*obj.L;
            m_26 = 0.5*A(1,1)*obj.L^2;
            m_33 = A(1,1)*obj.L;
            m_35 = -0.5*A(2,2)*obj.L^2;
            m_44 = A(3,3)*obj.L;
            m_46 = 0.5*A(1,3)*obj.L^2;
            m_55 = (1/3)*A(2,2)*obj.L^3;
            m_66 = (1/3)*A(1,1)*obj.L^3;

            M_over_rho = [0 0    0    0    0    0
                          0 m_22 m_23 m_24 0    m_26
                          0 0    m_33 0    m_35 0
                          0 0    0    m_44 0    m_46
                          0 0    0    0    m_55 0
                          0 0    0    0    0    m_66 ];
        end
        
        function theta_c = ThetaContact(obj)
            %ThetaContact: returns tilt angle (rad) for probe top to
            %contact the side wall of its melt channel.
            %
            %   Author: Nathan Daniel
            %   Date: 06/26/2020
            
            switch(obj.cross)
                case 'circle'
                    R = obj.r + obj.R_delta;
                    psi = acos(R/sqrt(obj.r^2 + obj.L^2));
                case 'square'
                    R = (obj.a + obj.R_delta)*[1 %First entry is distance to edge
                                               sqrt(2)]; %Second entry is distance to corner 
                    psi = acos(R/sqrt(obj.a^2 + obj.L^2));
                otherwise
                    error("%s is not a valid value for property 'cross'\n...(must be 'circle' or 'square')",obj.cross);
            end 
            theta_c = pi/2 - obj.gamma - psi;
                    
            %Negative angles will only occur if R<r (should never happen)
            %or R=r (trivial, no tilt possible)
            theta_c = max(theta_c,zeros(size(theta_c)));

            %Angles greater than pi/2 will only occur once the probe has passed
            %horizontal,at which point this analysis has not yet been assessed to
            %work (6-26-20).
            theta_c = min(theta_c,ones(size(theta_c))*pi/2);
        end
        
        function [A,B,lambda] = WaterPocketPitchDynamics_Linear(obj,T)
            %WaterPocketPitchDynamics_Linear provides the A matrix for the
            %pitch dynamics of a long, thin melt probe while it is fully
            %embedded in the liquid pocket that forms around it.
            %
            %   The input T (temperature in degrees C) is used to
            %   calculate the density (rho) and dynamic viscosity (mu) from
            %   arrays that are read in from a file. The gravitational
            %   acceleration is also read in from the file.
            %   
            %   TODO: account for Beta factor (vacuum contact with melt
            %   surface)
            %
            %   A modified inertia is calculated using the Added Mass
            %   tensor (see the AddedMassTensor method of this calss).
            %
            %   An A matrix is calculated using the following
            %   linearizations about theta=0:
            %       - cos(theta) (moment arm calculation) --> 1
            %           TODO: no skin friction term, is it even needed?
            %       - theta_d^2 (in skin friction term) --> 0
            %   These linearizations are very valid since theta will be
            %   very small (unlikely to be greater than 5 degrees)
            %
            %   A B matrix is also given which maps allows an input torque
            %   by dividing by the modified inertia value.
            %
            %   These dynamics are ONLY valid when theta is less than the
            %   theta required to contact the side walls (calculated by the
            %   ThetaContact method of this class), at which point static
            %   contact forces and melting dynamics become critical).
            
            %TODO: This should be fluid & planet independent
            load('meltProbeData.mat','mu_wat','rho_wat','T_wat','g_e');
            g = g_e;
            mu = interp1(T_wat,mu_wat,T);
            rho = interp1(T_wat,rho_wat,T);

            %I_yy = obj.I(2,2) + rho*obj.M_rho(5,5); %Added mass not valid
            %since fluid is heavily bounded
            I_yy = obj.I(2,2);
            
            A = [0                                            1
                 g*(obj.m*obj.CG - rho*obj.V*obj.CB)/I_yy    -obj.B1_mu*mu/I_yy ];
            B = [0
                 1/I_yy ];
            lambda = eig(A);
        end
        
        function [y,x_d] = WaterPocketPitchDynamics_Nonlinear(obj,x,u,T)
            load('meltProbeData.mat','mu_wat','rho_wat','T_wat','g_e');
            g = g_e;
            mu = interp1(T_wat,mu_wat,T);
            rho = interp1(T_wat,rho_wat,T);
            
            %I_yy = obj.I(2,2) + rho*obj.M_rho(5,5); %Added mass not valid
            %since fluid is heavily bounded
            I_yy = obj.I(2,2);
            
            x_d = [x(2);
                   (obj.m*g)*(obj.CG*sin(x(1)))/I_yy + ... %Gravity moment about drill tip
                   -(rho*g*pi*(obj.r^2)*obj.L)*(obj.CB*sin(x(1)))/I_yy + ... %Buoyancy moment about drill tip
                   -(obj.B1_mu*mu/I_yy)*x(2) + ... %Viscous damping moment about drill tip
                   -sign(x(2))*(obj.B2_rho*rho/I_yy)*x(2)^2 + ... %Pressure drag/skin friction moment about drill tip
                   u/I_yy]; %Disturbing torque
            y = [];
        end
        
        function M = IntegratedPressureDragMoment(obj,rho,nu,theta_d,approx_level)
            %Models Cd of cylinder with dependence on Reynolds number,
            %with Cd(Re) = { A*Re^b,     Re<1e2
            %                Cd100,      1e2<=Re<1e5 }
            % (Approx. 1: taken from http://ftp.demec.ufpr.br/disciplinas/TM045/FUNDAMENTALS_OF_AERODYNAMICS.PDF))
            %Expressing this as a function of the x position on the
            %vehicle,
            %with Cd(x) = { A*(c*theta_d*x/nu)^b,     x < x2
            %               Cd100,                      x2 < x < x5 }
            %where x2 = 1e2*nu/(c*theta_d) and x5 = 1e5*nu/(c*theta_d) and
            %c is the hydrodynamic chord length (c = diameter for a
            %cylinder).
            %Using this function to integrate along the pressure drag along
            %the length of the probe yields:
            % TODO
            
            %Cylindrical constants (taken from http://ftp.demec.ufpr.br/disciplinas/TM045/FUNDAMENTALS_OF_AERODYNAMICS.PDF)
            c = obj.d;
            A = 12.4286; b = -0.507; Cd100 = 1.2; %This is just an approximate fit
            
            x2 = 1e2*nu/(c*abs(theta_d)); %x at which Re = 1e2
            x5 = 1e5*nu/(c*abs(theta_d)); %x at which Re = 1e5
            
            switch(approx_level)
                case 1
                    if obj.L > x5 %Approx 1 not valid outside this range
                        error('Model for Cd(Re) for cylinder was not verified for Re>1e5, which has occured in this code.');
                    elseif obj.L > x2
                        M = 0.5*rho*(theta_d^2)*obj.r*(A*((obj.r*abs(theta_d)/nu)^b)*((1/(b+4))*(x2^(b+4)) + ((obj.r^2)/(b+2))*(x2^(b+2))) + ...
                                                       Cd100*(0.25*(obj.L^4 - x2^4) + ((obj.r^2)/2)*(obj.L^2 - x2^2)));
                    else
                        M = 0.5*rho*(theta_d^2)*obj.r*A*((obj.r*abs(theta_d)/nu)^b)*((1/(b+4))*(obj.L^(b+4)) + ((obj.r^2)/(b+2))*(obj.L^(b+2))); 
                    end
                case 2 %Since r<<L, assume it is 0 for calculating moment arm (r^2*L^2 terms --> 0 compared to L^4)
                    if obj.L > x5 %Approx 1 not valid outside this range
                        error('Model for Cd(Re) for cylinder was not verified for Re>1e5, which has occured in this code.');
                    elseif obj.L > x2
                        M = 0.5*rho*(theta_d^2)*obj.r*(A*((obj.r*abs(theta_d)/nu)^b)*(1/(b+4))*(x2^(b+4)) + ...
                                                       Cd100*(0.25*(obj.L^4 - x2^4)));
                    else
                        M = 0.5*rho*(theta_d^2)*obj.r*A*((obj.r*abs(theta_d)/nu)^b)*(1/(b+4))*(obj.L^(b+4));
                    end
            end
            M = -abs(M)*sign(theta_d);
        end
        
        function [P_req,P_dn,P_lat,L_crit] = RequiredPower(obj,ROP,T_amb,method)
            V_mps = ROP/3600;
            
            T_f = (7.4+11.4/2); %Observation in RECAS testing (i.e. condition dependent, but dependencies are unknown)
            
            load('meltProbeData.mat','g_e','T_ice','T_wat','rho_ice','rho_wat','cp_ice','cp_wat','k_ice','k_wat','mu_wat','L_melt','T_melt');
            rho_i = interp1(T_ice,rho_ice,(T_amb + T_melt)/2); %Suggestion from Ulamec: use averages over temp range from T_ice to T_melt
            rho_w = interp1(T_wat,rho_wat,(T_f + T_melt)/2);
            cp_i = interp1(T_ice,cp_ice,(T_amb + T_melt)/2);
            cp_w = interp1(T_wat,cp_wat,(T_f + T_melt)/2);
            k_i = interp1(T_ice,k_ice,(T_melt + T_amb)/2);
            k_w = interp1(T_wat,k_wat,(T_f + T_melt)/2);
            mu = interp1(T_wat,mu_wat,(T_f + T_melt)/2);
            alpha_i = k_i./(cp_i.*rho_i);
            alpha_w = k_w./(cp_w.*rho_w);
            
            switch method
                case 'Ulamec' %Ulamec 2006 paper, pure melt probe
                    alpha = 932; beta = 0.726;
                    
                    %Model only valid for certain values of x
                    x = obj.L./(V_mps*obj.A_x/pi);
%                     if (min(x) < 5e4)
%                         error('Max velocity is too large for model to be valid');
%                     end
%                     if (max(x) > 1e8)
%                         error('Min velocity is too small for model to be valid');
%                     end
                    P_dn = V_mps.*rho_i.*obj.A_x.*(L_melt + cp_i.*(T_melt - T_amb));
                    P_lat = (obj.A_x/pi)*V_mps.*(T_melt - T_amb)*alpha.*x.^beta;
                    P_req = P_dn + P_lat;
                    
                    P_dn((x<5e4)|(x>1e8)) = NaN;
                    P_lat((x<5e4)|(x>1e8)) = NaN;
                    P_req((x<5e4)|(x>1e8)) = NaN;
                    L_crit = obj.L*ones(size(P_req)); %This calculation estimates the lateral power required
                                    %such that the critical refreezing
                                    %length is the length of the probe.
                    
                case 'Talalay' %Talalay 2014 paper, Li 2020 paper on RECAS
                    k1 = 1;%1.04; %1.03 - 1.05
                    k2 = 1;%0.75; %0.7 - 0.8
                    k3 = 1;%1.03; %1.02 - 1.05 %TODO Not sure what this should be
                    
                    %Specific Energy Consumption (similar to cutting
                    %energy, but relative to descent amount not volume
                    Q12 = k1*rho_i.*obj.A_x.*(L_melt + k3*cp_i.*(T_melt - T_amb));
                    Q3 = k1*rho_i.*obj.A_x.*cp_w.*(T_f - T_melt);
                    
                    %Total power
                    P_req = V_mps.*(Q12+Q3)/k2;
                    P_dn = V_mps.*Q12/k2;
                    P_lat = V_mps.*Q3/k2;
                    L_crit = zeros(size(P_req)); %TODO: Double check what this calc should be
                    
                case 'Kowalski' %Kowalsky 2018 paper
                    alpha = 932; beta = 0.726; %called n & d in paper, but alpha & beta in Ulamec
                    
                    h_m_star = L_melt + cp_i.*(T_melt - T_amb); %reduced latent heat of melting
                    Q_c = obj.A_x*rho_i.*V_mps.*h_m_star; %Heat flowing from water film to ice
                    
                    %TODO: This term SHOULD also account for things like
                    %IceMole's ice screw (see paper)
                    F_star = (obj.m - rho_w*obj.V)*g_e; %buoyancy-corrected force
                    gamma_kow = (1./(20*alpha_w)).*(((rho_i./rho_w).*V_mps.*obj.r).^(4/3)).*((3*pi*mu./(2*F_star)).^(1/3));
                    
                    Q_h = ((7*gamma_kow+1)./(1-3*gamma_kow)).*Q_c;
                    Q_e = Q_h - Q_c;
                    
                    P_dn = Q_c;
                    P_req = Q_h;
                    P_lat = Q_h - Q_e; %Via convective losses
                    L_crit = ((Q_c.*(V_mps.^(beta-1)).*(obj.r^(2*(beta-1)))./(alpha*(T_melt - T_amb))).*((7*gamma_kow+1)./(1-3*gamma_kow) - 1)).^(1/beta);
                
                case 'Li' %Li 2020 paper
                    Q_i = obj.A_x*rho_i.*V_mps.*(L_melt + cp_i.*(T_melt - T_amb)); %Power to melt ice
                    gamma_li = 0.9; %0.5 to 0.9 (assumption) %Efficiency of thermal head
                    Q_t = Q_i/gamma_li; %Input power to thermal head
                    Q_w = Q_t - Q_i; %Power loss
                    
                    %If anyone ever reads this, I just want to say how
                    %absolutely goddamn stunned I am that Newton's method
                    %works perfectly for this ridiculous complex
                    %calculation comprised of multiple integrals of
                    %nonelementary functions. I mean math is crazy.
                    
                    %Use Newton's method to solve for critical refreezing length
                    if numel(V_mps)==1 %(single point or only T_amb is varied)
                        [L_crit,~,~] = newton(@(L) (Q_w - obj.QL_AamotAndLi(L,V_mps,T_amb)), obj.L*ones(size(T_amb)));
                    else %(Mesh or only V_mps is varied)
                        [L_crit,~,~] = newton(@(L) (Q_w - obj.QL_AamotAndLi(L,V_mps,T_amb)), obj.L*ones(size(V_mps))); 
                    end
                        
                        
                    P_dn = Q_i;
                    if L_crit>obj.L %Wasted power is enough to keep from refreezing
                        P_req = Q_t;
                        P_lat = Q_w;
                    else
                        P_req = Q_t + 0;
                        P_lat = Q_w + 0;
                    end
            end
        end
        
        function Q_L_star = QL_AamotAndLi(obj,L_star,V_mps,T_amb)
            
            %Ice properties
            load('meltProbeData.mat','T_ice','k_ice','rho_ice','cp_ice','T_melt');
            
            rho_i = interp1(T_ice,rho_ice,(T_amb + T_melt)/2); %Suggestion from Ulamec: use averages over temp range from T_ice to T_melt
            cp_i = interp1(T_ice,cp_ice,(T_amb + T_melt)/2);
            k_i = interp1(T_ice,k_ice,(T_melt + T_amb)/2);
            alpha_i = k_i./(cp_i.*rho_i);
            
            %Constant outside integral
            %(Note: The T_melt-T_amb is not mentioned in Li, but IS mentioned in Aamot)
            K = (8*k_i.*(T_amb-T_melt)/pi);
            err = 100/norm(K); %For err = n/abs(K), accurate to within n Watts (10 is the max accuracy I have been able to get to converge)
            
            %This integrand is described both in Li (2020) and Aamot(1967) (both
            %originate from Aamot which originates from analysis of cylindrical
            %heat flux)
            integrand = @(x,z,alpha,V) exp((-alpha.*(log(x).^2).*z)./V)./...
                                       (-x.*(-log(x)).*(besselj(0,obj.r.*(-log(x))).^2 + bessely(0,obj.r.*(-log(x))).^2));
            
            if (numel(V_mps)==1)&&(numel(T_amb)==1)
                integrand_i = @(x,z) integrand(x,z,alpha_i,V_mps);
                Q_L_star = K.*integral2(integrand_i,1e-10,1-1e-10,0,L_star,'AbsTol',err);
            elseif (numel(V_mps)>1)&&(numel(T_amb)==1)
                Q_L_star = zeros(size(V_mps));
                for i=1:numel(V_mps)
                    integrand_i = @(x,z) integrand(x,z,alpha_i,V_mps(i));
                    Q_L_star(i) = K.*integral2(integrand_i,1e-10,1-1e-10,0,L_star(i),'AbsTol',err);
                end
            elseif (numel(V_mps)==1)&&(numel(T_amb)>1)
                Q_L_star = zeros(size(T_amb));
                for i=1:numel(T_amb)
                    integrand_i = @(x,z) integrand(x,z,alpha_i(i),V_mps);
                    Q_L_star(i) = K(i).*integral2(integrand_i,1e-10,1-1e-10,0,L_star(i),'AbsTol',err);
                end
             elseif (numel(V_mps)==numel(T_amb))
                fprintf('Beginnning next newton iteration\n')
                L_star
                Q_L_star = zeros(size(V_mps));
                for i=1:numel(V_mps)
                    fprintf('\tCalculating L* for ROP = %f m/hr, T_amb = %f\n',V_mps(i)*3600,T_amb(i))
                    integrand_i = @(x,z) integrand(x,z,alpha_i(i),V_mps(i));
                    Q_L_star(i) = K(i).*integral2(integrand_i,1e-10,1-1e-10,0,L_star(i),'AbsTol',err);
                end
            else
                error('V_mps & T_amb have incorrect dimensions.')
            end
        end
    end
end



     

