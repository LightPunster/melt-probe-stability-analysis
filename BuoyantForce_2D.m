function F_b = BuoyantForce_2D(probe,rho,g,varargin)
%BuoyantForce_2D calculates the buoyant force on a melt probe.
%   Provided the 'probe' structure, which should contain
%       -The probe radius, probe.r
%       -The probe length, probe.L
%
%   And the additional parameters:
%       -rho = density of liquid (water)
%       -g = gravitational acceleration
%
%   Calculates the buoyant force experienced by the probe in the vertical
%   frame.
%
%   This may be calculated using Archimedes' Principle (AP), but
%   integrating the pressure loads on the probe yields an additional term
%   which is dependent on the quality of contact between the drill head and
%   the ice below the probe. This term produces no net torque, but does
%   reduces the upward buoyant force and create an additional horizontal
%   force in the direction that the probe is pointing.
%
%   In order to account for this force, input an additional three arguments:
%       - Beta = quality of contact factor (0<=Beta<=1)
%       - D = depth of probe melt head/drill in fluid
%           (NOTE/TODO: In the future, this may be replaced by a simple
%           pressure input to account for the fact that the refreezing ice
%           above the probe has a lower density than water).
%       -theta_phi = theta + phi angle of the probe w.r.t. the gravity vector

    switch nargin
        case 3
            F_b = rho*g*pi*(probe.r^2)*probe.L*[0;1];
        case 6
            Beta = varargin{1};
            D = varargin{2};
            theta_phi = varargin{3};
            F_b = rho*g*pi*(probe.r^2)*(probe.L*[0;1] + Beta*D*[sin(theta_phi);-cos(theta_phi)]);
        otherwise
            error('Expected 3 or 6 arguments, got %d.',nargin);
    end
end