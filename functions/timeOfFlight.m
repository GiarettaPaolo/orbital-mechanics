function [dt] = timeOfFlight(a, e, th_i, th_f, mu)

% timeOfFlight.m -
%
% PROTOTYPE:
% [dt] = timeOfFlight(a, e, th_i, th_f, mu)
%
% DESCRIPTION:
% Time elapsed between two points on the same orbit, as a function of the
% respective true anomalies
%
% INPUT:
% a                         [1x1]           Semi-major axis                         [km]
% e                         [1x1]           Eccentricity                            [-]
% th_i                      [1x1]           Initial true anomaly                    [rad]
% th_f                      [1x1]           Final true anomaly                      [rad]
% mu                        [1x1]           Gravitational parameter                 [km^3/s^2]
% OUTPUT:
% dt                        [1x1]           Time elapsed                            [km]

% Normalizing angles in [0, 2*pi) interval
th_i = mod(th_i, 2 * pi);
th_f = mod(th_f, 2 * pi);

if nargin == 4
    mu = 398600.433;
end

% Check for elliptical orbit
if (e >= 1) || ( e < 0)
    error('Orbit is not elliptical')
end

% 1. Eccentric anomaly from true anomaly
cosE_i = (e + cos(th_i)) / (1 + e*cos(th_i));
sinE_i = sqrt(1 - e^2) * sin(th_i) / (1 + e*cos(th_i));

cosE_f = (e + cos(th_f)) / (1 + e*cos(th_f));
sinE_f = sqrt(1 - e^2) * sin(th_f) / (1 + e*cos(th_f));

E_i = mod(atan2(sinE_i, cosE_i), 2 * pi);
E_f = mod(atan2(sinE_f, cosE_f), 2 * pi);

% 2. Mean anomaly from eccentric anomaly (Kepler's Equation)
dM = (E_f - E_i) - e * (sinE_f - sinE_i);

%3. Time elapsed from mean anomaly
dt = dM * sqrt(a^3 / mu);

% Check wether the final anomaly is smaller than the initial so a full
% period is added to avoid negative time intervals. 
if th_f < th_i
    dt = dt + 2*pi*sqrt(a^3 / mu);
end


end
