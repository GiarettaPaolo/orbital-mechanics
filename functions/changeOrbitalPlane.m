function [dv, om_f, th_f] = changeOrbitalPlane(a_i, e_i, i_i, OM_i, om_i, i_f, OM_f, mu, flag)

% changeOrbitalPlane.m - Change orbital plane of orbit.
%
% PROTOTYPE:
% [dv, om_f, th_f] = changeOrbitalPlane(a_i, e_i, i_i, OM_i, om_i, i_f, OM_f, mu, flag)
%
% DESCRIPTION:
% Changes the orbital plane of initial orbit to a
% designated orbital plane through a one-impulse maneuver. The shape of the
% orbit remains unchanged. The two possible maneuver points could result in
% a different delta-v cost; the flag argument can be used to change
% between those.
%
% INPUT:
% a_i                         [1x1]           Semi-major axis                         [km]
% e_i                         [1x1]           Eccentricity                            [-]
% i_i                         [1x1]           Initial inclination                     [rad]
% OM_i                        [1x1]           Initial RAAN                            [rad]
% om_i                        [1x1]           Initial anomaly of pericentre           [rad]
% i_f                         [1x1]           Final inclination                       [rad]
% OM_f                        [1x1]           Final RAAN                              [rad]
% mu                          [1x1]           Gravitational parameter                 [km^3/s^2]
% flag                        [1x1]           Maneuver point flag                     [-]
% OUTPUT:
% dv                          [1x1]           Burn delta-v                            [km/s]
% om_f                        [1x1]           Final anomaly of pericentre             [rad]
% th_f                        [1x1]           True anomaly of maneuver point          [rad]    

if (e_i > 1) || (e_i < 0)
    error('Orbit is not elliptical')
end

if nargin == 7
    flag = 1;
    mu = 398600.433;
elseif nargin == 8
    flag = 1;
end

% 1. Calculation of rotation angle alfa
dOM = OM_f - OM_i;
di = i_f - i_i;
alfa = acos( cos(i_i) * cos(i_f) + sin(i_i) * sin(i_f)*cos(dOM));

% 2. Four cases as a function of signs of dOM, di

if (dOM >= 0) && (di >= 0)
    
    cosu_i = (cos(alfa)*cos(i_i) - cos(i_f)) / sin(alfa) / sin(i_i);
    sinu_i = sin(i_f)*sin(dOM)/sin(alfa);
    
    cosu_f = (cos(i_i)-cos(alfa)*cos(i_f)) / sin(alfa)/sin(i_f);
    sinu_f = sin(i_i) * sin(dOM) / sin(alfa);
    
elseif (dOM >= 0) && (di <= 0)
    
    cosu_i = -(cos(alfa)*cos(i_i) - cos(i_f)) / sin(alfa) / sin(i_i);
    sinu_i = -sin(i_f)*sin(dOM)/sin(alfa);
    
    cosu_f = -(cos(i_i)-cos(alfa)*cos(i_f)) / sin(alfa)/sin(i_f);
    sinu_f = -sin(i_i) * sin(dOM) / sin(alfa);
    
elseif (dOM < 0) && (di >= 0)

cosu_i = (cos(alfa)*cos(i_i) - cos(i_f)) / sin(alfa) / sin(i_i);
sinu_i = sin(i_f)*sin(dOM)/sin(alfa);

cosu_f = (cos(i_i)-cos(alfa)*cos(i_f)) / sin(alfa)/sin(i_f);
sinu_f = sin(i_i) * sin(dOM) / sin(alfa);
    
elseif (dOM < 0) && (di <= 0)

cosu_i = -(cos(alfa)*cos(i_i) - cos(i_f)) / sin(alfa) / sin(i_i);
sinu_i = -sin(i_f)*sin(dOM)/sin(alfa);

cosu_f = -(cos(i_i)-cos(alfa)*cos(i_f)) / sin(alfa)/sin(i_f);
sinu_f = -sin(i_i) * sin(dOM) / sin(alfa);
    
end

% 3. Determination of true latitudes, arguments of pericentre and true
% anomaly
u_i = mod(atan2(sinu_i, cosu_i), 2 * pi);
u_f = mod(atan2(sinu_f, cosu_f), 2 * pi);
    
th_f = mod(u_i - om_i, 2 * pi);
om_f = u_f - th_f;
om_f = mod(om_f, 2 * pi);

% Manouver points on opposite side of the intersection line
% Final orbit parameters do not change but the delta-v cost can be
% different as it is only a function of trasverse velocity.
if flag == 0
    th_f = mod(th_f + pi, 2 * pi);  
end

% 4. Delta-v calculation dv = 2 * sqrt(mu / p) * v_th
dv = 2 * sin(alfa/2) * sqrt(mu / a_i / (1 - e_i^2)) * (1 + e_i * cos(th_f));

end




    
    