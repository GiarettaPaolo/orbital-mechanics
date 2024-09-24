function [a_t, e_t, th_t2] = tangent_orbit_equal_apocenter(ra_i, om_i, a_f, e_f, om_f, mu)

% tangent_orbit_equal_apocenter.m - Finds transfer orbit tangent at
% apocentre to the initial orbit and tangent to the second orbit
% orbit
%
% PROTOTYPE:
% [a_t, e_t, th_t2] = tangent_orbit_equal_apocenter(ra_i, om_i, a_f, e_f, om_f, mu)
%
% DESCRIPTION:
% Calculates the transfer orbit that is tangent at apocentre to the initial
% orbit and is tangent at another point to the final orbit. The problem can
% be expressed as a 4x4 system of non-linear equations in the eccentricity of
% transfer orbit, true anomaly of intersection point on final and transfer 
% orbits and transfer semi-major axis. The equations can be combined in a 
% 1x1 non linear equation that is solved by a fixed-point algorithm.

%
% INPUT:
% ra_i                 [1x1]          Initial apocentre                                 [km]
% om_i                 [1x1]          Initial anomaly of pericentre                     [rad]
% a_f                  [1x1]          Final semi-major axis                             [km]
% e_f                  [1x1]          Final eccentricity                                [-]
% om_f                 [1x1]          Final anomaly of pericentre                       [rad]
% mu                   [1x1]          Gravitational parameter                           [km^3/s^2]
%
% OUTPUT:
% a_t                  [1x1]          Transfer orbit semi-major axis                    [km]
% e_t                  [1x1]          Transfer orbit eccentricity                       [-]
% th_t2                [1x1]          True anomaly of intersection wrt transfer orbit   [km/s]

if nargin == 5
    mu = 398600.433;
end

% The equations are obtained by setting equal:
% the distances from the attractor at the tangency point
% p_t / (1 + e_t * cos(th_t2)) = p2 / (1 + e2 * cos(th2))
% The apocentre of the transfer orbit to the one given as input  
% ra_i = p_t / (1 + e_t)
% The flight path angles at the tangency point
% e_t * sin(th_t2) / (1 + e_t * cos(th_t2)) = e2 * sin(th2) / (1 + e2 * cos(th_t2))
% The true latitudes at the tangency point
% th_t2 + om_t = om2 + th2
% 
% The system reduces to the non linear equation K =(1 - e(th_t)) / e(th_t) * sin(th_t - dom) / sin(th_t)
% wehere K, e(th_t) are defined below

K = a_f * (1 - e_f^2) / ra_i / e_f;
dom = om_f - om_i;

e = @(th_t) e_f * sin(th_t - dom) / (sin(th_t) + e_f * sin(dom));
f = @(th_t) th_t + 0.05 * ( K -  (1 - e(th_t)) / e(th_t) * sin(th_t - dom) / sin(th_t));

% Fixed point iteration applied to function f to solve f(th) = th
% The weight 0.05 is to guarantee convergence 
[sol, ~] = ptofis(pi/2, f, 1000, 1e-12);
th_t2 = mod(sol, 2 * pi);

% Calculation of orbital parameters from converged solution
e_t = e(th_t2);
a_t = ra_i / (1 + e_t);

end

