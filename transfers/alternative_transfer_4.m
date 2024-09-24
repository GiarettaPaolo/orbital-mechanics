% Alternative transfer N°4
% The transfer consists of the following maneuvers:
% 1. On initial orbit wait to maneuver point.
% 2. Tangent burn at a suitable point.
% 3. Impulsive maneuver at apocentre such that transfer orbit is coplanar
% and tangent to the final orbit.
% 4. Impulsive burn at tangent intersection.
% 5. On final orbit wait to end point.
% Total dv = 2.4166 km/s
% Total time = 20587 s
% Remarks:
% i. The point on the initial orbit where the tangent burn is performed is
% chosen such that the resulting orbit has apocentre on the intersection
% betwen the initial and final orbital planes.
% ii. The semi-major axis of the second transfer obit has been left as a
% parameter a_t2. The optimal a_t2 has been found in order to minimize
% the total delta-v cost.
% iii. The final transfer orbit has been calculated by using the function 
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

% First bielliptic transfer orbit
a_t2 = 12235.0;      % 12235.0  
[e_t2, th_b1, dv_b1] = tangent_transfer_final_om(a1, e1, om1, om1+dom, a_t2, mu);
om_t2 = om1+dom;
i_t2 = i1;
OM_t2 = OM1;
th_t2 = th_b1 + om1 - om_t2;
dt_1 = timeOfFlight(a1, e1, th1, th_b1, mu);


% Plane change maneuver
[~, om_t3, th_b3] = changeOrbitalPlane(a_t2, e_t2, i_t2, OM_t2, om_t2, i2, OM2, mu, 1);
a_t3 = a_t2;
e_t3 = e_t2;
i_t3 = i2;
OM_t3 = OM2;
dt_2 = timeOfFlight(a_t2, e_t2, th_t2, th_b3);

% Tangent orbit 
[a_t4, e_t4, th_t4] = tangent_orbit_equal_apocenter(a_t3*(1 + e_t3), om_t3, a2, e2, om2, mu);
i_t4 = i_t3;
OM_t4 = OM_t3;
om_t4 = om_t3;
th_b4 = th_b3;
dt_3 = timeOfFlight(a_t4, e_t4, pi, th_t4, mu);

%Delta-v calculation
[~, v_b1] = kep2car(a_t2, e_t2, i_t2, OM_t2, om_t2, pi, mu);
[~, v_b2] = kep2car(a_t4, e_t4, i_t4, OM_t4, om_t4, pi, mu);
dv_b3 = norm(v_b2 - v_b1);

%Time of flight to get to final point
dt_4 = timeOfFlight(a2, e2, th_t4 + om_t4 - om2, th2, mu);

r4 = a_t4 * (1 - e_t4^2) / (1 + e_t4 * cos(th_t4));
dv_b4 = sqrt(2*mu*(1/r4 - 1/2/a2)) - sqrt(2*mu*(1/r4 - 1/2/a_t4));

dv_tot = abs(dv_b1) + abs(dv_b3) + abs(dv_b4)
T = dt_1 + dt_2 + dt_3 + dt_4

%% Plot Orbits

% Colors 
c1 = "#48E8C8";
c2 = "#FFBB00";
c3 = "#FF8100";

figure('Name','Alternative Transfer N4', 'NumberTitle', 'Off') 
set(gcf,'color','w');
Terra3d;

[X_t2, Y_t2, Z_t2] = plotOrbit([a_t2, e_t2, i_t2, OM_t2, om_t2, th_t2], mu, 2*pi, pi/1000);
[X_t4, Y_t4, Z_t4] = plotOrbit([a_t4, e_t4, i_t4, OM_t4, om_t4, th_t4], mu, 2 * pi, pi/1000);
plot3(X_t2, Y_t2, Z_t2, '--', 'Color', c1, 'LineWidth', 0.8)
plot3(X_t4, Y_t4, Z_t4,'--', 'Color', c2, 'LineWidth', 0.8)

[X1, Y1, Z1] = plotOrbit([a1, e1, i1, OM1, om1, th1], mu, 2*pi, pi/1000);
[X2, Y2, Z2] = plotOrbit([a2, e2, i2, OM2, om2, th2], mu, 2*pi, pi/1000);
[X_t2, Y_t2, Z_t2] = plotOrbit([a_t2, e_t2, i_t2, OM_t2, om_t2, th_t2], mu, th_b3 - th_t2, pi/1000);
[X_t3, Y_t3, Z_t3] = plotOrbit([a_t3, e_t3, i_t3, OM_t3, om_t3, th_b3], mu, 2*pi, pi/1000);
[X_t4, Y_t4, Z_t4] = plotOrbit([a_t4, e_t4, i_t4, OM_t4, om_t4, pi], mu, th_t4 - pi, pi/1000);

[aX_1, aY_1, aZ_1] = plotOrbit([a1, e1, i1, OM1, om1, 0], mu, pi, pi);
[aX_2, aY_2, aZ_2] = plotOrbit([a2, e2, i2, OM2, om2, 0], mu, pi, pi);
[aX_t2, aY_t2, aZ_t2] = plotOrbit([a_t2, e_t2, i_t2, OM_t2, om_t2, 0], mu, pi, pi);
[aX_t4, aY_t4, aZ_t4] = plotOrbit([a_t4, e_t4, i_t4, OM_t4, om_t4, 0], mu, pi, pi);

plot3(X1, Y1, Z1,'Color','b', 'LineWidth', 1)
plot3(X2, Y2, Z2,'Color','r', 'LineWidth', 1)
plot3(X_t2, Y_t2, Z_t2,'Color', c1, 'LineWidth', 1)
plot3(X_t4, Y_t4, Z_t4,'Color', c2, 'LineWidth', 1)

plot3(X1(1), Y1(1), Z1(1),'o','Color', 'b','MarkerSize', 6, 'MarkerFaceColor', 'b')
plot3(X2(1), Y2(1), Z2(1),'o','Color', 'r','MarkerSize', 6, 'MarkerFaceColor', 'r')
plot3(X_t2(1), Y_t2(1), Z_t2(1),'x','Color', c1,'MarkerSize', 8, 'LineWidth', 2)
plot3(X_t4(1), Y_t4(1), Z_t4(1),'x','Color', c2,'MarkerSize', 8, 'LineWidth', 2)
plot3(X_t4(end), Y_t4(end), Z_t4(end), 'x','Color', 'r','MarkerSize', 8, 'LineWidth', 2)

plot3(aX_1, aY_1, aZ_1, '-.','Color', 'b','LineWidth', 0.5)
plot3(aX_2, aY_2, aZ_2, '-.','Color', 'r','LineWidth', 0.5)
plot3(aX_t2, aY_t2, aZ_t2, '-.','Color', c1,'LineWidth', 0.5)
plot3(aX_t4, aY_t4, aZ_t4, '-.','Color', c2,'LineWidth', 0.5)

grid on
axis equal
box off
% title('\fontsize{15}{0}\selectfont Alternative Transfer N4','Interpreter','latex')
% subtitle_text = sprintf('dv = %.3f km/s, T = %0.f s', dv_tot, T);
% subtitle(subtitle_text, 'Interpreter', 'latex')
xlabel('x [km] - $\gamma$', 'Interpreter', 'latex') 
ylabel('y [km]', 'Interpreter', 'latex') 
zlabel('z [km]', 'Interpreter', 'latex') 
ax = gca;
ax.XRuler.Exponent = 0;
ax.YRuler.Exponent = 0;
ax.ZRuler.Exponent = 0;
