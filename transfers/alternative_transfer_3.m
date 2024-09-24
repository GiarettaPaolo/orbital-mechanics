% Alternative transfer N°3
% The transfer consists of the following maneuvers:
% 1. On initial orbit wait to maneuver point.
% 2. Change argument of pericentre maneuver to line up the absidal line with
% the intersection of the initial and final orbital planes.
% 3. Tangent burn at pericentre to desired transfer orbit.
% 4. Impulsive maneuver at apocentre such that transfer orbit is coplanar
% and tangent to the final orbit.
% 5. Impulsive burn at tangent intersection.
% 6. On final orbit wait to end point.
% Total dv = 2.5244 km/s
% Total time = 20175 s
% Remarks:
% i. The semi-major axis of the second transfer obit has been left as a
% parameter a_t2. The optimal a_t2 has been found in order to minimize
% the total delta-v cost.
% ii. The final transfer orbit has been calculated by using the function 
% tangent_orbit_equal_apocentre.m  
clear all

mu = 398600.433;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Initial and final orbits parameters %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r1 = [2254.3254; -8092.3126; -4199.8027];
v1 = [5.6120; 2.4220; -1.7020];

a2 = 16410.0000;
e2 = 0.2678;
i2 = 0.5612;
OM2 =0.4075;
om2 = 1.0700;
th2 = 1.3420;

[a1, e1, i1, OM1, om1, th1] = car2kep(r1, v1, mu);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Transfer algorithm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find the argument of perigee for transfer orbit
[~, ~, dom] = changeOrbitalPlane(a1, e1, i1, OM1, om1, i2, OM2, mu, 0);

% Change argument of pericentre
[dv_b1, om_t1, th_t1] = changePeriapsisArg(a1, e1, om1, dom, mu, 0);
a_t1 = a1;
e_t1 = e1;
i_t1 = i1;
OM_t1 = OM1;
dt_1 = timeOfFlight(a_t1, e_t1, th1, dom/2 + pi, mu);

% Second transfer orbit
a_t2 = 12223.8;         % <- optimal value ~12223.8
e_t2 = 1 - a_t1 * (1 - e_t1) / a_t2;
i_t2 = i_t1;
OM_t2 = OM_t1;
om_t2 = om_t1;
th_t2 = 0;
dv_b2 = sqrt(mu)*(sqrt((1 + e_t2) / (1 - e_t2) / a_t2) - sqrt((1 + e_t1) / (1 - e_t1) / a_t1));
dt_2 = timeOfFlight(a_t1, e_t1, -dom/2 + pi, 0, mu);
dt_3 = pi * sqrt( a_t2 ^ 3 / mu);

% Plane change maneuver
[~, om_t3, th_b3] = changeOrbitalPlane(a_t2, e_t2, i_t2, OM_t2, om_t2, i2, OM2, mu, 1);
a_t3 = a_t2;
e_t3 = e_t2;
i_t3 = i2;
OM_t3 = OM2;

% Tangent orbit 
[a_t4, e_t4, th_t4] = tangent_orbit_equal_apocenter(a_t3*(1 + e_t3), om_t3, a2, e2, om2, mu);
i_t4 = i_t3;
OM_t4 = OM_t3;
om_t4 = om_t3;
th_b4 = th_b3;
dt_4 = timeOfFlight(a_t4, e_t4, pi, th_t4, mu);
dt_5 = timeOfFlight(a2, e2, th_t4 + om_t4 - om2, th2, mu);

%Delta-v calculation
[~, v_b1] = kep2car(a_t2, e_t2, i_t2, OM_t2, om_t2, pi, mu);
[~, v_b2] = kep2car(a_t4, e_t4, i_t4, OM_t4, om_t4, pi, mu);
dv_b3 = norm(v_b2 - v_b1);

r4 = a_t4 * (1 - e_t4^2) / (1 + e_t4 * cos(th_t4));
dv_b4 = sqrt(2*mu*(1/r4 - 1/2/a2)) - sqrt(2*mu*(1/r4 - 1/2/a_t4));

dv_tot = abs(dv_b1) + abs(dv_b2) + abs(dv_b3) + abs(dv_b4)
T = dt_1 + dt_2 + dt_3 + dt_4 + dt_5

%% Plot Orbits

% Colors 
c1 = "#48E8C8";
c2 = '#42AA79';
c3 = "#FF8100";
c4 = "#FFBB00";

figure('Name','Alternative Transfer N3', 'NumberTitle', 'Off') 
set(gcf,'color','w');
Terra3d;

[X_t1, Y_t1, Z_t1] = plotOrbit([a_t1, e_t1, i_t1, OM_t1, om_t1, th_t1], mu, 2*pi, pi/1000);
[X_t2, Y_t2, Z_t2] = plotOrbit([a_t2, e_t2, i_t2, OM_t2, om_t2, th_t2], mu, 2 * pi, pi/1000);
% [X_t3, Y_t3, Z_t3] = plotOrbit([a_t3, e_t3, i_t3, OM_t3, om_t3, 0], mu, 2*pi , pi/1000);
[X_t4, Y_t4, Z_t4] = plotOrbit([a_t4, e_t4, i_t4, OM_t4, om_t4, 0], mu, 2*pi , pi/1000);
plot3(X_t1, Y_t1, Z_t1, '--', 'Color', c1, 'LineWidth', 0.8)
plot3(X_t2, Y_t2, Z_t2,'--', 'Color', c2, 'LineWidth', 0.8)
% plot3(X_t3, Y_t3, Z_t3,'--', 'Color', c3, 'LineWidth', 0.8)
plot3(X_t4, Y_t4, Z_t4,'--', 'Color', c4, 'LineWidth', 0.8)

[X1, Y1, Z1] = plotOrbit([a1, e1, i1, OM1, om1, th1], mu, 2*pi, pi/1000);
[X2, Y2, Z2] = plotOrbit([a2, e2, i2, OM2, om2, th2], mu, 2*pi, pi/1000);
[X_t1, Y_t1, Z_t1] = plotOrbit([a_t1, e_t1, i_t1, OM_t1, om_t1, th_t1], mu, 2*pi - th_t1, pi/1000);
[X_t2, Y_t2, Z_t2] = plotOrbit([a_t2, e_t2, i_t2, OM_t2, om_t2, 0], mu, pi, pi/1000);
[X_t3, Y_t3, Z_t3] = plotOrbit([a_t3, e_t3, i_t3, OM_t3, om_t3, th_b3], mu, 2*pi, pi/1000);
[X_t4, Y_t4, Z_t4] = plotOrbit([a_t4, e_t4, i_t4, OM_t4, om_t4, pi], mu, th_t4 - pi, pi/1000);

[aX_1, aY_1, aZ_1] = plotOrbit([a1, e1, i1, OM1, om1, 0], mu, pi, pi);
[aX_2, aY_2, aZ_2] = plotOrbit([a2, e2, i2, OM2, om2, 0], mu, pi, pi);
[aX_t1, aY_t1, aZ_t1] = plotOrbit([a_t1, e_t1, i_t1, OM_t1, om_t1, 0], mu, pi, pi);
[aX_t2, aY_t2, aZ_t2] = plotOrbit([a_t2, e_t2, i_t2, OM_t2, om_t2, 0], mu, pi, pi);
[aX_t4, aY_t4, aZ_t4] = plotOrbit([a_t4, e_t4, i_t4, OM_t4, om_t4, 0], mu, pi, pi);

plot3(X1, Y1, Z1,'Color','b', 'LineWidth', 1)
plot3(X2, Y2, Z2,'Color','r', 'LineWidth', 1)
plot3(X_t1, Y_t1, Z_t1,'Color', c1, 'LineWidth', 1)
plot3(X_t2, Y_t2, Z_t2,'Color', c2, 'LineWidth', 1)
plot3(X_t4, Y_t4, Z_t4,'Color', c4, 'LineWidth', 1)

plot3(X1(1), Y1(1), Z1(1),'o','Color', 'b','MarkerSize', 6, 'MarkerFaceColor', 'b')
plot3(X2(1), Y2(1), Z2(1),'o','Color', 'r','MarkerSize', 6, 'MarkerFaceColor', 'r')
plot3(X_t1(1), Y_t1(1), Z_t1(1),'x','Color', c1,'MarkerSize', 8, 'LineWidth', 2)
plot3(X_t2(1), Y_t2(1), Z_t2(1),'x','Color', c2,'MarkerSize', 8, 'LineWidth', 2)
plot3(X_t4(1), Y_t4(1), Z_t4(1),'x','Color', c3,'MarkerSize', 8, 'LineWidth', 2)
plot3(X_t4(end), Y_t4(end), Z_t4(end), 'x','Color', 'r','MarkerSize', 8, 'LineWidth', 2)

plot3(aX_1, aY_1, aZ_1, '-.','Color', 'b','LineWidth', 0.5)
plot3(aX_2, aY_2, aZ_2, '-.','Color', 'r','LineWidth', 0.5)
plot3(aX_t1, aY_t1, aZ_t1, '-.','Color', c1,'LineWidth', 0.5)
plot3(aX_t2, aY_t2, aZ_t2, '-.','Color', c2,'LineWidth', 0.5)
plot3(aX_t4, aY_t4, aZ_t4, '-.','Color', c4,'LineWidth', 0.5)

grid on
axis equal
box off
% title('\fontsize{15}{0}\selectfont Alternative Transfer N3','Interpreter','latex')
% subtitle_text = sprintf('dv = %.3f km/s, T = %0.f s', dv_tot, T);
% subtitle(subtitle_text, 'Interpreter', 'latex')
xlabel('x [km] - $\gamma$', 'Interpreter', 'latex') 
ylabel('y [km]', 'Interpreter', 'latex') 
zlabel('z [km]', 'Interpreter', 'latex') 
ax = gca;
ax.XRuler.Exponent = 0;
ax.YRuler.Exponent = 0;
ax.ZRuler.Exponent = 0;
