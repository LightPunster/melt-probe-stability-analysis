function [B1_mu,B2_rho] = RotationalDragCoeffs(probe,varargin)
%RotationalDrag calculates coefficients of viscous and pressure drag for
%rotation of a melt probe. Normalized by density to allow later dependence
%on temperature.
%   Viscous Drag: Tau_d1 = B1_over_rho*rho*omega;
%   Pressure Drag/Skin Friction: Tau_d2 = B2_over_rho*rho*omega^2


%method = 'Heiss & Coull (1952)';
method = 'Bowen & Masliyah (1973)';

switch probe.cross
    case 'circle'
        %For C_d = the coefficient of drag of a cross section exposed to
        %lateral flow, the torque may be found by integrating the drag
        %equation along the length of the probe
        %dTau = 0.5*rho*V(l)^2*Cd*l*dA (differential torque on one dl-long slice
        %Tau = integral{(0.5*rho*V(l)^2)*Cd*l}dA
        %    = integral(0,L){(0.5*rho*(l^2)*(theta_d^2)*Cd*l*r}dl
        %    = 0.5 * rho * theta_d^2 * r * integral(0,L){(l^3)*Cd}dl
        %Some papers give Cd=1.2 at high Reynolds numbers (>=10e4)
        %   (http://ftp.demec.ufpr.br/disciplinas/TM045/FUNDAMENTALS_OF_AERODYNAMICS.PDF)
        %   (http://sv.20file.org/up1/916_0.pdf)
        %    = 0.125 * Cd * r * L^4
        %But lower Reynolds numbers show much higher values approximately:
        %   (http://ftp.demec.ufpr.br/disciplinas/TM045/FUNDAMENTALS_OF_AERODYNAMICS.PDF)
        %   Cd = 12.4286*Re^-0.507 for Re<10^2, and ~1.2 above that
        %      = 12.4286*(V*d/nu)^-0.507
        %      = 12.4286 * (d/nu)^-0.507 * (theta_d*l)^-0.507
        %      = 12.4286 * (d/nu)^-0.507 * (theta_d^-0.507) * l^-0.507
        %Making the full integral:
        %Tau = 0.5 * 12.4286 * rho * nu^0.507 * theta_d^1.493 * r * (2*r)^-0.507 * integral(0,L1){l^2.493}dl +
        %      0.125 * 1.2 * r * (L^4 - L1^4)
        %Where L1 is the L at which R = 100 ==> L1 = 100*nu/(theta_d*2*r)
        %    = 6.2143 * rho*nu^0.507 * theta_d^1.493 * 0.7037 * r^0.493 * (1/3.493) * (100*nu/theta_d*2*r)^3.493 + 
        %      0.125 * 1.2 * r * (L^4 - (100*nu/theta_d*2*r)^4)
        %    = 1.0767e+06 * rho* nu^4 * theta_d^-2 * r^-3 + 0.15*r*L^4 -
        %    937500*theta_d^-4
        %      0.125 * 1.2 * r * (L^4 - (100*nu/theta_d*2*r)^4)
        %Getting VERY messy: TODO
        
        %This has actually shown to be EXTREMELY accurate for VERNE, but
        %need to test for other shapes & at high speeds, or at least arrive
        %at by approximation & prove error. HOWEVER, that is a good "TODO"
        %for later or even outside the scope of this project.
        B2_rho = 0.125*DragCoeff(probe.cross,'Anderson')*probe.r*probe.L^4; %Pressure drag effects

        %Calculated a factors w.r.t. Stoke's Law drag on sphere (https://www.tandfonline.com/doi/pdf/10.1080/02786828708959128
        switch(method)
            case 'Heiss & Coull (1952)'
                d_v = 2*((3/4)*probe.L.*probe.r.^2).^(1/3); %Diameter of sphere with same volume
                d_n = 2*probe.r; %Diameter of sphere with same projected area normal to dir of motion
                Psi = (4*pi*(d_v/2).^2) / (2*pi*probe.r.*probe.L + 2*pi*probe.r.^2); %Surface area of sphere with same volume/surface area of probe
                K = 10^(-0.25*((Psi*d_v/d_n)^0.5)*((d_v/d_n)-1) + log10((d_v/d_n)*Psi^0.5)); %Factor to Stoke's Law drag
                
                B1_mu = 3*pi*d_v*K*(probe.L/2); %L/2 is very rough estimate for moment arm
            case 'Bowen & Masliyah (1973)'
                d_p = (2*probe.L + 4*probe.r)./pi; %Diameter of sphere with same perimeter
                S = (2*pi*probe.r.*probe.L + 2*pi*probe.r.^2) / (4*pi*d_p.^2); %Suface area of probe / surface area of perimeter-equivalent sphere
                K_p = 0.392 + 0.621*S - 0.040*S^2;
                
                B1_mu = 3*pi*d_p*K_p*(probe.L/2); %L/2 is very rough estimate for moment arm
          
            
            otherwise
                error('No valid method for calculating viscous drag coefficient of crossection "%s" provided.',probe.cross)
        end
    case 'square'
        B2_rho = 0;
        B1_mu = 0;

    otherwise
        error('No valid cross section provided')
end

end

function Cd = DragCoeff(cross,source)
    switch(cross)
        case 'circle'
            switch(source)
                case 'Anderson' %http://ftp.demec.ufpr.br/disciplinas/TM045/FUNDAMENTALS_OF_AERODYNAMICS.PDF
                    Cd = 1.2; %@R = 10^5
                case 'Brushi' %http://sv.20file.org/up1/916_0.pdf
                    Cd = 1.16; %average of results for R ~(10^4, 10^5)
                case 'min'
                    Cd = min([1.2,1.16]);
                case 'max'
                    Cd = max([1.2,1.16]);
                case 'avg'
                    Cd = mean([1.2,1.16]);
            end
    end
end