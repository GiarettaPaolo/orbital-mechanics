function [X, Y, Z] = plotOrbit( kepEl , mu, deltaTh , stepTh)

% plotOrbit.m - 3D coordinates of orbital points for plotting.
%
% PROTOTYPE:
% [X, Y, Z] = plotOrbit( kepEl , mu, deltaTh , stepTh)
%
% DESCRIPTION:
% Returns the coordinates with respect to the ECI frame of a orbit. The
% points are parametrized by true anomaly. The spanned interval and
% increment can be chosen. 
%
% INPUT:
% kepEl                    [1x6]          Keplerian parameters              
% mu                       [1x1]          Gravitational parameter           [km^3/s^2]
% deltaTh                  [1x1]          True anomaly range                [rad]
% stepTh                   [1x1]          True anomaly step increment       [rad]
%
% OUTPUT:
% X                        [1xn]          X coordinate of points            [km]   
% Y                        [1xn]          Y coordinate of points            [km]  
% Z                        [1xn]          Z coordinate of points            [km] 

% Set default values for deltaTh, stepTh [deltaTh = 2pi, stepTh = pi/180]
if nargin == 2
    deltaTh = 2*pi;
    stepTh = pi / 180;
elseif nargin == 3
    stepTh = pi / 180;
end


if deltaTh < 0
    deltaTh = mod(deltaTh, 2 * pi);
end
    
a = kepEl(1);
e = kepEl(2);
i = kepEl(3);
OM = kepEl(4);
om = kepEl(5);
th_start = kepEl(6);

% True anomalies vector
th_vect = [0: stepTh: deltaTh] + th_start;
n = length(th_vect);

% Pre-allocating for speed
X = zeros(1, n);
Y = zeros(1, n);
Z = zeros(1, n);

% Get cartesian coordinates of plot points
for k=1:n
    
    [r, ~] = kep2car(a, e, i, OM, om, th_vect(k), mu);
    X(k) = r(1);
    Y(k) = r(2);
    Z(k) = r(3);
    
end








