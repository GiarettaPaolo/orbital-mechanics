function [e_f, th_b1, om_t, dv] = tangent_transfer_coaxial_orbits(kepEl1, kepEl2, a_f, mu)

% tangent_transfer_coaxial_orbits - Change orbit shape while remaining tangent
% to initial orbit.
% 
%
% PROTOTYPE:
% [e_f, th_b1, om_t, dv_b1] = tangent_transfer_coaxial_orbits(kepEl1, kepEl2, a_f, mu)
%
% DESCRIPTION:
% Performs an impulsive burn tangent to the initial orbit. The
% final semi-major axis is given. The algorithm finds the burn point for which
% the orbit, obtained after a subsequent plane change maneuver, is
% coplanar and coaxial with the final orbit.
%
% INPUT:
% kepEl1               [1x6]          Keplerian parameters of initial point
% kepEl2               [1x6]          Keplerian parameters of final point              
% a_f                  [1x1]          Final semi-major axis                     [km]
% mu                   [1x1]          Gravitational parameter                   [km^3/s^2]
%
% OUTPUT:
% e_f                  [1x1]          Final eccentricity                        [-]
% th_b1                [1x1]          True anomaly of burn point                [rad]
% om_t                 [1x1]          Anomaly of pericentre of transfer orbit   [rad]
% dv                   [1x1]          Delta-v of impulsive burn                 [km/s]


if nargin == 3
    mu = 398600.433;
end

a1 = kepEl1(1);
e1 = kepEl1(2);
i1 = kepEl1(3);
OM1 = kepEl1(4);
om1 = kepEl1(5);

i2 = kepEl2(3);
OM2 = kepEl2(4);
om2 = kepEl2(5);

% Fixed-point iteration algorithm to find the true anomaly of maneuver
% point to have desired final om after the plane change maneuver.
% The error is estimated through succesive differences.

th_b1 = 0;          % Initial guess of anomaly of burn point
nmax = 1000;        % Maximum number of iterations
tol = 1e-12;        % Desired tolerance
it = 0;             % Iteration counter
err = tol + 1;
while (err > tol) && (it < nmax)
    
    [e_f, om_t, dv] = tangent_burn(a1, e1, om1, th_b1, a_f, mu);
    [~, om_f, ~] = changeOrbitalPlane(a_f, e_f, i1, OM1, om_t, i2, OM2, mu, 1); 
    
    % Fixed point iteration, weight of 0.1 to guarantee convergence
    th_b1 = th_b1 - 0.1 * (om_f - om2);
    err = abs(om_f - om2);
    it = it + 1;
    
end
th_b1 = mod( th_b1, 2 * pi);

end