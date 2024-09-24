function [dv_b1, dv_b2, th_f, dt] = bitangent_changeOrbitShape(a_i, e_i, om_i, a_f, e_f, om_f, mu, flag)

% bitangent_changeOrbitShape.m - Bitangent transfer between coaxial
% elliptical orbits.
%
% PROTOTYPE:
% [dv_b1, dv_b2, th_f, dt] = bitangent_changeOrbitShape(a_i, e_i, om_i, a_f, e_f, om_f, mu, flag)
%
% DESCRIPTION:
% Bitangent transfer orbit between two coaxial elliptical orbits. The
% transfer is a two-impulse maneuver such that the burns are performed at
% absidal points. Depending on the relative position of initial and
% transfer orbits four cases are possible(pericentre - pericentre,
% apocentre - apocentre, pericentre - apocentre, apocentre - pericentre).
% The flag argument is 1 if the first impulse is carried at the pericentre of
% initial orbit, and 0 for apocentre.
%
% INPUT:
% a_i                      [1x1]          Initial semi-major axis           [km]
% e_i                      [1x1]          Initial eccentricity              [-]
% om_i                     [1x1]          Initial anomaly of pericentre     [rad]
% a_f                      [1x1]          Final semi-major axis             [km]
% e_f                      [1x1]          Final eccentricity                [-]
% om_f                     [1x1]          Final anomaly of pericentre       [rad]
% mu                       [1x1]          Gravitational parameter           [km^3/s^2]
% flag                     [1x1]          Flag argument                     [-]
%
% OUTPUT:
% dv_b1                    [1x1]          First impulse delta-v             [km/s]
% dv_b2                    [1x1]          Second impulse delta-v            [km/s]
% th_f                     [1x1]          Final true anomaly                [rad]

% Setting mu and flag default values
if nargin == 6
    flag = 1;
    mu = 398600.344;
elseif nargin == 7
    flag = 1;
end

% Calculating apocentre and pericentre distances
rp_i = a_i * (1 - e_i);
ra_i = a_i * (1 + e_i);
rp_f = a_f * (1 - e_f);
ra_f = a_f * (1 + e_f);

if om_i == om_f
    
    if flag == 1
        
        % From perigee to apogee
        a12 = (rp_i + ra_f ) / 2;
        th_f = pi;
        dv_b1 = sqrt(mu * ( 2 / rp_i - 1 / a12 )) - sqrt(mu * ( 2 / rp_i - 1 / a_i )) ;
        dv_b2 =  sqrt(mu * ( 2 / ra_f - 1 / a_f )) - sqrt(mu * ( 2 / ra_f - 1 / a12 ));
        dt = pi * sqrt( a12 ^ 3 / mu);
        
    elseif flag == 0
        
        % From apogee to perigee
        a12 = (ra_i + rp_f ) / 2;
        th_f = pi;
        dv_b1 = sqrt(mu * ( 2 / ra_i - 1 / a12 )) - sqrt(mu * ( 2 / ra_i - 1 / a_i ));
        dv_b2 = sqrt(mu * ( 2 / rp_f - 1 / a_f )) - sqrt(mu * ( 2 / rp_f - 1 / a12 ));
        dt = pi * sqrt( a12 ^ 3 / mu);
    end
    
elseif abs(om_i - om_f) == pi
    
    if flag == 1
        
        % From perigee to apogee
        a12 = (rp_i + rp_f ) / 2;
        th_f = pi;
        dv_b1 = sqrt(mu * ( 2 / rp_i - 1 / a12 )) - sqrt(mu * ( 2 / rp_i - 1 / a_i ));
        dv_b2 =  sqrt(mu * ( 2 / rp_f - 1 / a_f )) - sqrt(mu * ( 2 / rp_f - 1 / a12 ));
        dt = pi * sqrt( a12 ^ 3 / mu);
        
    elseif flag == 0
        
        % From apogee to perigee
        a12 = (ra_i + ra_f ) / 2;
        th_f = pi;
        dv_b1 = sqrt(mu * ( 2 / ra_i - 1 / a12 )) - sqrt(mu * ( 2 / ra_i - 1 / a_i )) ;
        dv_b2 = sqrt(mu * ( 2 / ra_f - 1 / a_f )) - sqrt(mu * ( 2 / ra_f - 1 / a12 ));
        dt = pi * sqrt( a12 ^ 3 / mu);
    end
    
else
    error('Orbits are not coaxial')
end
    
end

    
