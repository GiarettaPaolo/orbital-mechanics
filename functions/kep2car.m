function [r, v] = kep2car(a, e, i, OM, om, th, mu)

% kep2car.m - 
%
% PROTOTYPE:
% [r, v] = kep2car(a, e, i, OM, om, th, mu)
%
% DESCRIPTION:
% Conversion from Keplerian elements to Cartesian coordinates. Angles in
% radians.
%
% INPUT:
% a                       [1x1]           Semi-major axis                      [km]
% e                       [1x1]           Eccentricity                         [-]
% i                       [1x1]           Inclination                          [rad]
% OM                      [1x1]           RAAN                                 [rad]
% om                      [1x1]           Pericentre anomaly                   [rad]
% th                      [1x1]           True anomaly                         [rad]
% mu                      [1x1]           Gravitational parameter              [km^3/s^2]
% OUTPUT:
% r                       [3x1]           Position vector                      [km]
% v                       [3x1]           Velocity vector                      [km/s]
%

% 1. Derived orbital parameters
p = a * (1 - e ^ 2);
r_norm = p / ( 1 + e * cos(th) );

% 2. State vector in the perifocal frame 
r_pf = r_norm .* [cos(th); sin(th); 0];
v_pf = sqrt( mu / p) .* [-sin(th);( e + cos(th)); 0];

% 3. Rotation matrix from perifocal to ECI coordinates
R = zeros(3);
R(1, 1) = cos(om)*cos(OM) - sin(om)*cos(i)*sin(OM);
R(1, 2) = -sin(om)*cos(OM)-cos(om)*cos(i)*sin(OM);
R(1, 3) = sin(i)*sin(OM);
R(2, 1) = cos(om)*sin(OM) + sin(om)*cos(i)*cos(OM);
R(2, 2) = -sin(om)*sin(OM) + cos(om)*cos(i)*cos(OM);
R(2, 3) = -sin(i)*cos(OM);
R(3, 1) = sin(om)*sin(i);
R(3, 2) = cos(om)*sin(i);
R(3, 3) = cos(i);

% 4. position and velocity in ECI (change of basis)
r = R * r_pf;
v = R * v_pf;
