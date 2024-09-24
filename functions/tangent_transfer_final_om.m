function [e_f, th_b1, dv] = tangent_transfer_final_om(a_i, e_i, om_i, om_f, a_f, mu)

% tangent_transfer_final_om.m - Change orbit shape while remaining tangent
% to initial orbit.
% orbit
%
% PROTOTYPE:
% [e_f, th_b1, dv_b1] = tangent_transfer_final_om(a_i, e_i, om_i, om_f, a_f, mu)
%
% DESCRIPTION:
% Performs an impulsive burn tangent to the initial orbit. The
% final semi-major axis and argument of pericentre are given. The solution
% is obtained through a fixed-point iteration algorithm.
%
% INPUT:
% a_i                  [1x1]          Initial semi-major axis           [km]
% e_i                  [1x1]          Initial eccentricity              [-]
% om_i                 [1x1]          Initial anomaly of pericentre     [rad]
% om_f                 [1x1]          Final anomaly of pericentre       [rad]
% a_f                  [1x1]          Final semi-major axis             [km]
% mu                   [1x1]          Gravitational parameter           [km^3/s^2]
%
% OUTPUT:
% e_f                  [1x1]          Final eccentricity                [-]
% th_b1                [1x1]          True anomaly of burn point        [rad]
% dv                   [1x1]          Delta-v of impulsive burn         [km/s]

if nargin == 3
    mu = 398600.433;
end

% Fixed-point iteration algorithm to find the true anomaly of maneuver
% point to have desired final om. The error is estimated through succesive
% differences.

th_b1 = 0;          % Initial guess of anomaly of burn point
nmax = 1000;        % Maximum number of iterations
tol = 1e-12;        % Desired tolerance
it = 0;             % Iteration counter
err = tol + 1;
while (err > tol) && (it < nmax)
    
    [e_f, om_t, dv] = tangent_burn(a_i, e_i, om_i, th_b1, a_f, mu);   
    
    % Fixed point iteration, weight of 0.1 to guarantee convergence
    th_b1 = th_b1 - 0.1 * (om_t - om_f);
    err = abs(om_t - om_f);
    it = it + 1;
    
end
th_b1 = mod( th_b1, 2 * pi);

end