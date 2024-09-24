function [a, e, i OM, om, th] = car2kep(r, v, mu)

% car2kep.m - Conversion from Cartesian coordinates to Keplerian elements
%
% PROTOTYPE:
% [a, e, i, OM, om, th] = car2kep(r, v, mu)
%
% DESCRIPTION:
% Conversion from Cartesian coordinates to Keplerian elements. Angles in
% radians.
%
% INPUT:
% r                        [3x1]          Position vector                   [km]
% v                        [3x1]          Velocity vector                   [km/s]
% mu                       [1x1]          Gravitational parameter           [km^3/s^2]
%
% OUTPUT:
% a                        [1x1]          Semi-major axis                   [km]
% e                        [1x1]          Eccentricity                      [-]
% i                        [1x1]          Inclination                       [rad]
% OM                       [1x1]          RAAN                              [rad]
% om                       [1x1]          Pericentre anomaly                [rad]
% th                       [1x1]          True anomaly                      [rad]


% 1. Magnitude position and velocity vectors
r_norm = norm(r);                                 
v_norm = norm(v);                                

% 2. Specific angular momentum
h = cross(r, v);                                      
h_norm = norm(h);                               

% 3. Inclination
i = acos(h(3) / h_norm);                       

% 4. Eccentricity and eccentricity vector
vr = dot(r, v);  % radial speed
e_vect = ( (v_norm ^ 2 - mu / r_norm) * r - vr * v ) / mu;
e = norm(e_vect);

% 5. Semi-major axis
a = mu / (2 * mu / r_norm - v_norm ^ 2);

% 6. Node vector
N = cross([0;0;1], h);
N_norm = norm(N);

% 7. RAAN
if N(2) >= 0  % Check which branch [0, pi] - (pi, 2*pi) to choose
    OM = acos(N(1) / N_norm);
else
    OM = 2*pi - acos(N(1) / N_norm);
end
    
% 8. Pericenter anomaly
if e_vect(3) >= 0  % Check which branch [0, pi] - (pi, 2*pi) to choose
    om = acos(dot(N, e_vect) / (e * N_norm));
else
    om = 2*pi - acos(dot(N, e_vect) / (e * N_norm));
end
    
% 9. True anomaly
if vr >= 0  % Check which branch [0, pi] - (pi, 2*pi) to choose
    th = acos(dot(r, e_vect) / (e * r_norm));
else
    th = 2*pi - acos(dot(r, e_vect) / (e * r_norm));
end
    
    
