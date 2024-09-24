% extra transfer 2
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


% Plane Change manouver 
[dv_b1, om_t1, th_t1] = changeOrbitalPlane(a1, e1, i1, OM1, om1, i2, OM2, mu, 1);
a_t1 = a1;
e_t1 = e1;
OM_t1 = OM2;
i_t1 = i2;

% Time transfer to get to maneuver point 1
th_b1 = th_t1;
dt_1 = timeOfFlight(a1, e1, th1, th_b1);


% Change Periapsis Argument  maneuver
dom = om2 - om_t1;
[dv_b2, om_t2, th_t2] = changePeriapsisArg(a_t1, e_t1, om_t1, dom,mu, 1);
a_t2 = a_t1;
e_t2 = e_t1;
OM_t2 = OM_t1;
i_t2 = i_t1;

% Time transfer to get to maneuver point 2
th_b2 = dom / 2 + 2 * pi;
dt_2 = timeOfFlight(a_t1, e_t1, th_t1, th_b2);

% Time transfer to maneuver point 3 (perigee)
th_b3 = 0;
dt_3 = timeOfFlight(a_t2, e_t2, th_t2, th_b3);

% Bitangent transfer pericenter - apocenter
rb = 2 * a2 * (1 + e2);
[dv_b3, dv_b4, dv_b5, th_b5, dt_4] = bielliptic_changeOrbitShape(a_t2, e_t2, om_t2, a2, e2,om2, rb, mu, 1);
a_t3 = (a_t2 * (1 - e_t2) + rb) / 2;
e_t3 = (rb - a_t2 * (1 - e_t2)) / 2 / a_t3;
OM_t3 = OM2;
om_t3 = om2;
i_t3 = i2;

a_t4 = (a2 * (1 - e2) + rb) / 2;
e_t4 = (rb - a2 * (1 - e2)) / 2 / a_t4;
OM_t4 = OM2;
om_t4 = om2;
i_t4 = i2;

% Time transfer to get to final point
dt_5 = timeOfFlight(a2, e2, th_b5, th2, mu);

% Total time
T = dt_1 + dt_2 + dt_3 + dt_4 + dt_5;
% Total delta V
dv_tot = abs(dv_b1) + abs(dv_b2) + abs(dv_b3) + abs(dv_b4) + abs(dv_b5);


%% Plot Orbits

[X1, Y1, Z1] = plotOrbit([a1, e1, i1, OM1, om1, th1], mu, 2*pi, pi/100);
[X2, Y2, Z2] = plotOrbit([a2, e2, i2, OM2, om2, th2], mu, 2*pi, pi/100);
[X_t1, Y_t1, Z_t1] = plotOrbit([a_t1, e_t1, i_t1, OM_t1, om_t1, th_t1], mu, th_b2 - th_t1, pi/100);
[X_t2, Y_t2, Z_t2] = plotOrbit([a_t2, e_t2, i_t2, OM_t2, om_t2, th_t2], mu, th_b3 - th_t2, pi/100);
[X_t3, Y_t3, Z_t3] = plotOrbit([a_t3, e_t3, i_t3, OM_t3, om_t3, 0], mu, pi, pi/100);
[X_t4, Y_t4, Z_t4] = plotOrbit([a_t4, e_t4, i_t4, OM_t4, om_t4, pi], mu, pi, pi/100);

[aX_1, aY_1, aZ_1] = plotOrbit([a1, e1, i1, OM1, om1, 0], mu, pi, pi);
[aX_2, aY_2, aZ_2] = plotOrbit([a2, e2, i2, OM2, om2, 0], mu, pi, pi);
[aX_t1, aY_t1, aZ_t1] = plotOrbit([a_t1, e_t1, i_t1, OM_t1, om_t1, 0], mu, pi, pi);
[aX_t2, aY_t2, aZ_t2] = plotOrbit([a_t2, e_t2, i_t2, OM_t2, om_t2, 0], mu, pi, pi);
[aX_t3, aY_t3, aZ_t3] = plotOrbit([a_t3, e_t3, i_t3, OM_t3, om_t3, 0], mu, pi, pi);

Terra3d;
plot3(X1, Y1, Z1, 'b', 'LineWidth', 1)
plot3(X2, Y2, Z2, 'r', 'LineWidth', 1)
plot3(X_t1, Y_t1, Z_t1, 'c', 'LineWidth', 1)
plot3(X_t2, Y_t2, Z_t2, 'k', 'LineWidth', 1)
plot3(X_t3, Y_t3, Z_t3, 'g', 'LineWidth', 1)
plot3(X_t4, Y_t4, Z_t4, 'g', 'LineWidth', 1)

plot3(X1(1), Y1(1), Z1(1), 'ob', 'LineWidth', 2)
plot3(X2(1), Y2(1), Z2(1), 'or', 'LineWidth', 2)
plot3(X_t1(1), Y_t1(1), Z_t1(1), 'xc', 'LineWidth', 2)
plot3(X_t2(1), Y_t2(1), Z_t2(1), 'xk', 'LineWidth', 2)
plot3(X_t3(1), Y_t3(1), Z_t3(1), 'xg', 'LineWidth', 2)
plot3(X_t4(1), Y_t4(1), Z_t4(1), 'xg', 'LineWidth', 2)
plot3(X_t4(end), Y_t4(end), Z_t4(end), 'xr', 'LineWidth', 2)

plot3(aX_1, aY_1, aZ_1, '-.b', 'LineWidth', 0.5)
plot3(aX_2, aY_2, aZ_2, '-.r', 'LineWidth', 0.5)
plot3(aX_t1, aY_t1, aZ_t1, '-.c', 'LineWidth', 0.5)
plot3(aX_t3, aY_t3, aZ_t3, '-.g', 'LineWidth', 0.5)

grid on


