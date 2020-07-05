function R = Vert2Chan(theta_ch,varargin)
%Vert2Chan Generates rotation matrix for transforming from the vertical
%frame to the frame used to describe the melt channel
%   For the 2D version, the assumption is that all motion occurs in the
%   global x,z frame

switch nargin
    case 2
        dim = varargin{1};
    otherwise
        dim = 3;
end

switch(dim)
    case 2
        R = [ cos(pi/2+theta_ch)  sin(pi/2+theta_ch)
              -sin(pi/2+theta_ch) cos(pi/2+theta_ch) ];
    otherwise
        error('Invalid value for "dimension" parameter provided.');
end
end
