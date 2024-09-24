function [a_t, e_t, i_t, OM_t, om_t, th_t1, th_t2, dv_min] = direct_optimal_transfer(kepEl1, kepEl2, mu)

% direct_orbital_transfer - Computes the most efficient (delta-v wise)
% transfer orbit passing through initial and final points.
% 
%
% PROTOTYPE:
% [a_t, e_t, i_t, OM_t, om_t, th_t1, th_t2, dv_min] = direct_optimal_transfer(kepEl1, kepEl2, mu)
%
% DESCRIPTION:
% Given the initial and final point expressed by Keplerian parameters,
% the direct orbital transfer that requires the least delta-v is
% calculated. An infinitude of orbits pass through the assigned points and
% they can all be parametrized by one parameter (in this case the true
% anomaly of the first burn point). The algorithm returns the trajectory that 
% minimizes the delta-v required for the two impulsive maneuvers. A check
% is performed to prevent the satellite from intersecting the Earth's
% surface.
%
% INPUT:
% kepEl1               [1x6]          Keplerian parameters of initial point
% kepEl2               [1x6]          Keplerian parameters of final point              
% mu                   [1x1]          Gravitational parameter                        [km^3/s^2]
%
% OUTPUT:
% a_t                  [1x1]          Semi-major axis of transfer orbit              [-]
% e_t                  [1x1]          Transfer orbit eccentricity                    [-]
% i_t                  [1x1]          Inclination of transfer orbit                  [rad]
% OM_t                 [1x1]          RAAN of transfer orbit                         [rad]
% om_t                 [1x1]          Anomaly of pericentre of transfer orbit        [rad]
% th_t1                [1x1]          True anomaly of first point on transfer orbit  [rad]
% th_t2                [1x1]          True anomaly of second point on transfer orbit [rad]
% dv_min               [1x1]          Delta-v of optimal trajectory                  [km/s]

Re = 6378;

if nargin == 2
    mu = 398600.433;
end

a1 = kepEl1(1);
e1 = kepEl1(2);
i1 = kepEl1(3);
OM1 = kepEl1(4);
om1 = kepEl1(5);
th1 = kepEl1(6);

a2 = kepEl2(1);
e2 = kepEl2(2);
i2 = kepEl2(3);
OM2 = kepEl2(4);
om2 = kepEl2(5);
th2 = kepEl2(6);

% Default values if no solution is found
om_t = nan;
th_t1 = nan;
th_t2 = nan;
dv_min = inf;

[r1, v1] = kep2car(a1, e1, i1, OM1, om1, th1, mu);
[r2, v2] = kep2car(a2, e2, i2, OM2, om2, th2, mu);

% Two impulse secant approach

% Local coordinate axis
% x_loc is parallel do r1, z_loc to the specific angular momentum and y_loc
% constructed to have an orthonormal basis.
x_loc = r1 / norm(r1);
z_loc = cross(x_loc, r2);   

if z_loc(3) < 0          % <- <0 for prograde orbits >0 for retrograde orbits
    z_loc = -z_loc;
end

z_loc = z_loc / norm(z_loc);
y_loc = cross(z_loc, x_loc);
y_loc = y_loc / norm(y_loc);


% Determination of angle deltaTh between position vectors
if dot(r2, y_loc) >= 0
    deltaTh = acos( dot(r2, r1) / norm(r2) / norm(r1));
else
    deltaTh = 2 * pi - acos( dot(r2, r1) / norm(r2) / norm(r1));
end

% Inclination   
i_t = acos(z_loc(3));

% RAAN
N = cross([0;0;1], z_loc);
N = N / norm(N);
if N(2) >= 0
    OM_t = acos(N(1));
else
    OM_t = 2 * pi - acos(N(1));
end

% Main optimization loop (find orbit with minimum delta_v)
% The algorithm minimizes the total delta-v required by the transfer by
% parametrizing orbits as a function of the true anomaly of first point.
% At each iteration the interval [a, b] is divided into n sample points
% where the cost function is evalued at each point. The local minimum is
% calculated for the sample points and for the next iteration the interval
% [a, b] is centered at the minimum point. The algorithm halts when the 
% step between successive sample points is lower than the tolerance.

tol = 1e-12;              % Tolerance
n = 20;                   % Number of sample points for iteration
step = 1 + tol;             
a = 0;                    % Initial interval lower bound
b = 2*pi;                 % Initial interval lower bound
dv_min = inf;               
th_min = 0;               

while step > tol
    
    step = (b - a) / n;
    
    for th_t1 = a:step:b
        % Orbital parameters of transfer orbit
        th_t2 = mod(th_t1 +  deltaTh, 2 * pi);
        
        e_t = (norm(r1) - norm(r2)) / (norm(r2) * cos(th_t2) - norm(r1) * cos(th_t1));
        a_t = norm(r1) * (1 + e_t * cos(th_t1)) / ( 1 - e_t ^ 2);
        
        % Invalid orbits (either negative eccentricity or the satellite
        % would intersect the Earth's surface during transit)
        if e_t < 0
            continue
        elseif (a_t * (1 - e_t) < Re) && (th_t2 < th_t1)
            continue
        end
        

        % Argument of pericenter (from normalized eccentricity vector obtained 
        % by rotating the x_loc vector by -th_t1)
        e_unit = cos(th_t1) * x_loc - sin(th_t1) * y_loc;
        
        % Argument of pericentre of transfer orbit
        if e_unit(3) >= 0
            om_t = acos(dot(e_unit, N));
        else
            om_t = 2 * pi - acos(dot(e_unit, N));
        end
        
        % Delta-v calculated as norm of velocity vectors difference
        [~, v_t1] = kep2car(a_t, e_t, i_t, OM_t, om_t, th_t1, mu);
        [~, v_t2] = kep2car(a_t, e_t, i_t, OM_t, om_t, th_t2, mu);
        dv = norm(v_t2 - v2) + norm(v_t1 - v1);
        
        % Check for local optimal solution
        if dv < dv_min
            dv_min = dv;
            th_min = th_t1;
        end
        
    end
    
    % Updating the searching interval [a, b] for next iteration
    a = th_min - step;
    b = th_min + step;
    
end

end

