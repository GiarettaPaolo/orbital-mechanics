function [dv_1, dv_2, dv_3, th_f, dt] = bielliptic_changeOrbitShape(a_i, e_i, om_i, a_f, e_f, om_f, rb, mu, flag)

% bielliptic_changeOrbitShape.m - Bitangent transfer between coaxial
% elliptical orbits.
%
% PROTOTYPE:
% [dv_1, dv_2, dv_3, th_f, dt] = bielliptic_changeOrbitShape(a_i, e_i, om_i, a_f, e_f, om_f, rb, mu, flag)
%
% DESCRIPTION:
% Bielliptic transfer orbit between two coaxial elliptical orbits. The
% transfer is a three-impulse maneuver such that the burns are performed at
% absidal points. Depending on the relative position of initial and
% transfer orbits four cases are possible(pericentre - pericentre,
% apocentre - apocentre, pericentre - apocentre, apocentre - pericentre).
% The flag argument is 1 if the first impulse is carried at the pericentre of
% initial orbit, and 0 for apocentre. The parameter rb is the absidal
% distance on the first transfer orbit.
%
% INPUT:
% a_i                      [1x1]          Initial semi-major axis           [km]
% e_i                      [1x1]          Initial eccentricity              [-]
% om_i                     [1x1]          Initial anomaly of pericentre     [rad]
% a_f                      [1x1]          Final semi-major axis             [km]
% e_f                      [1x1]          Final eccentricity                [-]
% om_f                     [1x1]          Final anomaly of pericentre       [rad]
% rb                       [1x1]          Transfer absidal distance         [km]
% mu                       [1x1]          Gravitational parameter           [km^3/s^2]
% flag                     [1x1]          Flag argument                     [-]
%
% OUTPUT:
% dv_1                     [1x1]          First impulse delta-v             [km/s]
% dv_2                     [1x1]          Second impulse delta-v            [km/s]
% dv_3                     [1x1]          Third impulse delta-v             [km/s]
% th_f                     [1x1]          Final true anomaly                [rad]

if nargin == 6
    flag = 1;
    mu = 398600.433;
elseif nargin == 7
    flag = 1;
end

if om_i == om_f
    
    % pericentre - pericentre
    if flag == 1
        r1 = a_i * (1 - e_i);
        r3 = a_f * (1 - e_f);
        th_f = 0;
    % apocentre - apocentre    
    elseif flag == 0
        r1 = a_i * ( 1 + e_i);
        r3 = a_f * (1 + e_f);
        th_f = pi;
    
elseif abs(om_i - om_f) == pi
      
    % pericentre - apocentre
    if flag == 1
      r1 = a_i * (1 - e_i);
      r3 = a_f * (1 + e_f);
      th_f = pi;
    
    % apocentre - pericentre   
    elseif flag == 0
        r1 = a_i * ( 1 + e_i);
        r3 = a_f * (1 - e_f);
        th_f = pi;
    end
else
    error('orbite non coassiali')
    end

% delta-v values 
dv_1 = sqrt(2*mu*(1/r1 - 1/(r1+rb))) - sqrt(2*mu*(1/r1 - 1/2/a_i));
dv_2 = sqrt(2*mu*(1/rb - 1/(rb + r3))) - sqrt(2*mu*(1/rb - 1/(r1+rb)));
dv_3 = sqrt(2*mu*(1/r3 - 1/2/a_f)) - sqrt(2*mu*(1/r3 - 1/(r3 + rb)));
% Transfer period
dt = pi * sqrt( (r1 + rb) ^3 / 8 / mu) + pi * sqrt( (r3 + rb) ^3 / 8 / mu);
      
      
end