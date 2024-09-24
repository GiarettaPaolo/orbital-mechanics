% extra transfer 3
% Supplementary code / transfers that were abandoned 
clear all

mu = 398600.433;

r1 = [2254.3254; -8092.3126; -4199.8027];
v1 = [5.6120; 2.4220; -1.7020];

a2 = 16410.0000;
e2 = 0.2678;
i2 = 0.5612;
OM2 =0.4075;
om2 = 1.0700;
th2 = 1.3420;
[a1, e1, i1, OM1, om1, th1] = car2kep(r1, v1, mu);

% Find the argument of perigee for transfer orbit
[~, ~, dom] = changeOrbitalPlane(a1, e1, i1, OM1, om1, i2, OM2, mu, 0);

% Change argument of pericentre
[dv_b1, om_t1, th_t1] = changePeriapsisArg(a1, e1, om1, dom, mu, 0);
a_t1 = a1;
e_t1 = e1;
i_t1 = i1;
OM_t1 = OM1;
dt_1 = timeOfFlight(a_t1, e_t1, th1, dom/2 + pi, mu);

% First bielliptic transfer orbit
a_t2 = 20000;         
e_t2 = 1 - a_t1 * (1 + e_t1) / a_t2;
i_t2 = i_t1;
OM_t2 = OM_t1;
om_t2 = om_t1 + pi;
th_t2 = 0;
dv_b2 = sqrt(mu)*(sqrt((1 + e_t2) / (1 - e_t2) / a_t2) - sqrt((1 - e_t1) / (1 + e_t1) / a_t1));
dt_2 = timeOfFlight(a_t1, e_t1, dom/2 + pi, 0, mu);
dt_3 = pi * sqrt( a_t2 ^ 3 / mu);

% Plane change maneuver
[dv_b3, om_t3, th_b3] = changeOrbitalPlane(a_t2, e_t2, i_t2, OM_t2, om_t2, i2, OM2, mu, 0);
a_t3 = a_t2;
e_t3 = e_t2;
i_t3 = i2;
OM_t3 = OM2;

% Secant solution
p2 = a2 * (1 - e2^2);
p_t3 = a_t3 * (1 - e_t3^2);

%C = A * cos(th) + B * sin(th) = G cos(th + phi)
C = p2 - p_t3;
A = (e2 * p_t3 * cos(om2 - om_t3) - e_t3 * p2);
B = e2 * p_t3 * sin(om2 - om_t3);

G = sqrt(A^2 + B^2);
phi = atan2(-B/G, A/G);

th_b4 = mod(acos(C/G) - phi, 2 * pi); % <- put a - sign infront of the acos for the other (sub_optimal) secant solution
th_t4 = mod(th_b4 + om_t3 - om2, 2 * pi);

% Delta-v calculation
vr1 = sqrt(mu / p_t3) * e_t3 * sin(th_b4);
vt1 = sqrt(mu / p_t3) * (1 + e_t3 * cos(th_b4));
vr2 = sqrt(mu / p2) * e2 * sin(th_t4);
vt2 = sqrt(mu / p2) * (1 + e2 * cos(th_t4));
dv_b4 = sqrt((vr2 - vr1) ^ 2 + (vt2 - vt1) ^ 2);

dv_tot = abs(dv_b1) + abs(dv_b2) + abs(dv_b3) + abs(dv_b4)

%% Plot Orbits

[X1, Y1, Z1] = plotOrbit([a1, e1, i1, OM1, om1, th1], mu, 2*pi, pi/1000);
[X2, Y2, Z2] = plotOrbit([a2, e2, i2, OM2, om2, th2], mu, 2*pi, pi/1000);
[X_t1, Y_t1, Z_t1] = plotOrbit([a_t1, e_t1, i_t1, OM_t1, om_t1, th_t1], mu, 2 * pi, pi/1000);
[X_t2, Y_t2, Z_t2] = plotOrbit([a_t2, e_t2, i_t2, OM_t2, om_t2, th_t2], mu, th_b3 - th_t2, pi/1000);
[X_t3, Y_t3, Z_t3] = plotOrbit([a_t3, e_t3, i_t3, OM_t3, om_t3, th_b3], mu, th_b4 - th_b3, pi/1000);

[aX_1, aY_1, aZ_1] = plotOrbit([a1, e1, i1, OM1, om1, 0], mu, pi, pi);
[aX_2, aY_2, aZ_2] = plotOrbit([a2, e2, i2, OM2, om2, 0], mu, pi, pi);
[aX_t1, aY_t1, aZ_t1] = plotOrbit([a_t1, e_t1, i_t1, OM_t1, om_t1, 0], mu, pi, pi);
[aX_t2, aY_t2, aZ_t2] = plotOrbit([a_t2, e_t2, i_t2, OM_t2, om_t2, 0], mu, pi, pi);

Terra3d;
plot3(X1, Y1, Z1, 'b', 'LineWidth', 1)
plot3(X2, Y2, Z2, 'r', 'LineWidth', 1)
plot3(X_t1, Y_t1, Z_t1, 'c', 'LineWidth', 1)
plot3(X_t2, Y_t2, Z_t2, 'k', 'LineWidth', 1)
plot3(X_t3, Y_t3, Z_t3, 'g', 'LineWidth', 1)

plot3(X1(1), Y1(1), Z1(1), 'ob', 'LineWidth', 2)
plot3(X2(1), Y2(1), Z2(1), 'or', 'LineWidth', 2)
plot3(X_t1(1), Y_t1(1), Z_t1(1), 'xc', 'LineWidth', 2)
plot3(X_t2(1), Y_t2(1), Z_t2(1), 'xk', 'LineWidth', 2)
plot3(X_t3(1), Y_t3(1), Z_t3(1), 'xg', 'LineWidth', 2)
plot3(X_t3(end), Y_t3(end), Z_t3(end), 'xr', 'LineWidth', 2)

plot3(aX_1, aY_1, aZ_1, '-.b', 'LineWidth', 0.5)
plot3(aX_2, aY_2, aZ_2, '-.r', 'LineWidth', 0.5)
plot3(aX_t1, aY_t1, aZ_t1, '-.c', 'LineWidth', 0.5)

grid on


