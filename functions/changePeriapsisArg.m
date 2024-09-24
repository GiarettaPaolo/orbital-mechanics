function [dv, om_f, th_f] = changePeriapsisArg(a_i, e_i, om_i, dom, mu, flag)

%changePeriapsisArg.m - Rotates the initial orbit by changing the argument
%of pericentre.
%
% PROTOTYPE:
% [dv, om_f, th_f] = changePeriapsisArg(a_i, e_i, om_i, dom,mu, flag)
%
% DESCRIPTION:
% Rotates the orbit aorund an axis perpendicular to the orbital plane, by
% a one impulse maneuver. Both possible intersection points result in the
% same delta-v cost, the burn is carried at pericentre for flag=1 or
% apocentre for flag=0;
%
% INPUT:
% a_i                         [1x1]           Semi-major axis                         [km]
% e_i                         [1x1]           Eccentricity                            [-]
% om_i                        [1x1]           Initial anomaly of pericentre           [rad]
% dom                         [1x1]           Change in anomaly of pericentre         [rad]
% mu                          [1x1]           Gravitational parameter                 [km^3/s^2]
% flag                        [1x1]           Maneuver point flag                     [-]
% OUTPUT:
% dv                          [1x1]           Burn delta-v                            [km/s]
% om_f                        [1x1]           Final anomaly of pericentre             [rad]
% th_f                        [1x1]           True anomaly of maneuver point          [rad]       

if nargin == 4
    flag = 1;
    mu = 398600.433;
elseif nargin == 5
    flag = 1;
end

% delta-v calculated as 2 * sqrt(mu / p) * e * sin(dom/2)
dv = 2 * sqrt(mu / a_i / (1 - e_i ^ 2)) * e_i * sin(dom / 2);
om_f = mod(om_i + dom, 2 * pi);
% by changing the flag argument the other burn point is diametrically
% opposite. (change of pi of true anomaly)
th_f = mod(pi * (1 + flag) - dom/2, 2 * pi);


end

