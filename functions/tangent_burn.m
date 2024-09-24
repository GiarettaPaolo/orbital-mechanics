function [e_f, om_f, dv] = tangent_burn(a_i, e_i, om_i, th_i, a_f, mu)

% tangent_burn.m - Change orbit shape while remaining tangent to initial
% orbit
%
% PROTOTYPE:
% [e_f, om_f, dv] = tangent_burn(a_i, e_i, om_i, th_i, a_f, mu)
%
% DESCRIPTION:
% Performs a tangent burn at a specific point on the initial orbit. The
% final semi-major axis and true anomaly of maneuver point are given.
%
% INPUT:
% a_i                  [1x1]          Initial semi-major axis           [km]
% e_i                  [1x1]          Initial eccentricity              [-]
% om_i                 [1x1]          Initial anomaly of pericentre     [rad]
% th_i                 [1x1]          True anomaly of impulsive burn    [rad]
% a_f                  [1x1]          Final semi-major axis             [km]
% mu                   [1x1]          Gravitational parameter           [km^3/s^2]
%
% OUTPUT:
% e_f                  [1x1]          Final eccentricity                [-]
% om_f                 [1x1]          Final anomaly of pericentre       [rad]
% dv                   [1x1]          Burn delta-v                      [km/s]

if nargin == 5
    mu = 398600.433;
end

% Point on the initial orbit where the impulsive maneuver is performed
p_i = a_i * ( 1 - e_i ^ 2);
r =  p_i/ (1 + e_i * cos(th_i));

% Initial velocities (before impulsive maneuver)
vr_i = sqrt(mu / p_i) * e_i * sin(th_i);
vt_i = sqrt(mu / p_i) * ( 1 + e_i * cos(th_i));
v_i = sqrt(vr_i ^ 2 + vt_i ^ 2);

% Final velocities (after impulsive maneuver)
% The velocities are obtained scaling by v_f / v_i due to tangential burn
v_f = sqrt(2 * mu * (1 / r - 1 / 2 /a_f));
vt_f = vt_i * v_f / v_i;
vr_f = vr_i * v_f / v_i;
dv = v_f - v_i;

% Calculation of e_t from known semi-major axis and specific angular
% momentum ( p = sqrt(h * mu))
p_f = (r * vt_f) ^ 2 / mu;
e_f = sqrt(1 - p_f / a_f);

% Calculation of om_f from the final radial and tangential velocities
costh_f = (p_f / r - 1) / e_f;
sinth_f = vr_f / sqrt(mu / p_f) / e_f;

th_f = mod( atan2(sinth_f, costh_f), 2*pi);
om_f =mod( om_i + th_i - th_f, 2 * pi);

end