function [F_point,F_wall,status] = StaticContactForces_2D(probe,theta,phi,mu,F_in,varargin)
%StaticContactForces_2D calculates static contact forces for melt probe at
%wall and melt head/drill point
%   Provided the 'probe' object, which should contain
%       -The probe radius, probe.r
%       -The probe length, probe.L
%       -The probe diagonal, probe.diag = sqrt(probe.L^2 + probe.r^2)
%       -The probe diagonal angle, probe.gamma = atan(probe.r/probe.L)
%       -The diameter of the melt pocket, probe.R
%
%   And the additional parameters:
%       -theta = angle of the probe w.r.t. the melt channel
%       -phi = angle of the melt channel w.r.t. the gravity vector
%       -mu = coefficient of friction of melt channel walls
%       -F_in = sum of all other forces (gravity, buoyancy, etc.) on probe
%       as a 2D vector in the vertical frame at the CG
%       -Optional: T_in sum of additional torques NOT due to F_in
%
%   Returns the solution to the static force-moment balance problem as the 
%   vertical frame reaction forces at the point of contact on the side wall
%   (if it exists) and at the base.
%   Returns NaNs if no solution exists and an error if infinite solutions
%   exists. The solution state may be seen in the additional output
%   "status".

%Calculate theta contact and determine if contact is taking place
theta_c = probe.ThetaContact(); %Theta to contact side walls
theta = lim(theta,[-theta_c,theta_c]); %theta may not be greater than theta_c
contact = (theta==theta_c)||(theta==-theta_c); %If theta==theta_c, probe is contacting side walls

%Transform to probe frame
F_in_pf = [cos(theta+phi) sin(theta+phi); -sin(theta+phi) cos(theta+phi)]*F_in;

%Equations:
%   1. Moment balance in probe frame about drill point
%   2. Force balance along probe lateral axis
%   3. Force balance along probe longitudinal axis
%   4. Frictional relationship at side wall
A = [0 0        -cos(probe.gamma)*probe.diag*contact -sin(probe.gamma)*probe.diag*contact ;
     1 0         1                                    0;
     0 1         0                                    1;
     0 0         sin(theta)-mu*cos(theta)             cos(theta)+mu*sin(theta)];

switch nargin
    case 5
        b = [F_in_pf(1)*probe.CG; %Note: perpendicular component of force*CG gives all external torques due to forces
            -F_in_pf(1);
            -F_in_pf(2); 0];
    case 6
        T_in = varargin{1};
        b = [F_in_pf(1)*probe.CG + T_in;
            -F_in_pf(1);
            -F_in_pf(2); 0];
    otherwise
        error('Expected 5 or 6 arguments, got %d.',nargin);
end

[x,status,free,~] = linEqnSol(A,b);

%Get solution vector
if strcmp(status,'Unique solution')
    F_point_pf = [x(1); x(2)];
    F_wall_pf = [x(3); x(4)];
elseif strcmp(status,'No solution')
    F_point_pf = [NaN;NaN];
    F_wall_pf = [NaN;NaN];
else %infinite solutions, likely due to no contact with side wall (F_wall = 0)
    if sum(free==3) || sum(free==4) %If either component of the wall contact force is free, then there is no contact 
        x = x(:,1); %Then both components of the wall contact force are zero, and x can be set to the homogeneous solution
        F_point_pf = [x(1); x(2)];
        F_wall_pf = [x(3); x(4)];
    else
        error('Infinite solutions and unable to simplify; more constraints required.');
    end
end

%Transform back into vertical frame
F_point = [cos(theta+phi) -sin(theta+phi); sin(theta+phi) cos(theta+phi)]*F_point_pf;
F_wall = [cos(theta+phi) -sin(theta+phi); sin(theta+phi) cos(theta+phi)]*F_wall_pf;

end

