function [th_f] = anomalyChange(a, e, th_i, dt, mu, tol)

% anomalyChange.m - True anomaly after set time has elapsed.
%
% PROTOTYPE:
% [th_f] = anomalyChange(a, e, th_i, dt, mu, tol)
%
% DESCRIPTION:
% Returns the true anomaly of the satellite given the inital orbit, initial
% true anomaly and time elapsed.
%
% INPUT:
% a                         [1x1]           Semi-major axis                     [km]
% e                         [1x1]           Eccentricity                        [-]
% th_i                      [1x1]           Initial true anomaly                [rad]
% dt                        [1x1]           Time elapsed                        [s]
% mu                        [1x1]           Gravitational parameter             [km^3/s^2]
% OUTPUT:
% th_f                      [1x1]           Final true anomaly                  [km]


th_i = mod(th_i, 2 * pi);

if nargin == 4
    mu = 398600.433;
    tol = 1e-6;
elseif nargin == 5
    tol = 1e-6;
end

if (e >= 1) || ( e < 0)    
    error('Orbit is not elliptical')
end

% Initial time elapsed since periapsis (proportional to mean anomaly)
% Final time elapsed since periapsis (dt_i + dt)
dt_i  = timeOfFlight(a, e, 0, th_i, mu);
dt_f = dt_i + dt;

% Mean anomaly 
M = dt_f * sqrt(mu / a^3);

% Newton's method to solve Kepler's equation
E = M;              % Initial guess
it = 0;             % Iteration counter
nmax = 1000;        % Maximum number of iterations
while abs(E - e*sin(E) - M) > tol
    % E_n+1 = E_n - f(E_n) / f'(E_n)
    E = E - (E - e*sin(E) - M) / ( 1 -e * cos(E));
    it = it + 1;
    
    if it > nmax
        error('Newtons method does not converge')
    end
    
end

% True anomaly from eccentric anomaly
costh_f = (cos(E) - e) / (1 - e * cos(E));
sinth_f = sin(E) / sqrt(1 - e ^ 2) * (1 + e * costh_f);

% Normalizing converged result in the [0, 2 * pi) interval
th_f = mod(atan2(sinth_f, costh_f), 2 * pi);

end

