% Alternative transfer N°2
% The transfer consists of the following maneuvers:
% 1. On initial orbit wait to maneuver point.
% 2. Tangent burn at suitable point on initial orbit to get to desired
% transfer orbit.
% 3. Plane change maneuver at minimum delta-v cost point.
% 4. Transfer orbit and final orbits are tangent at apocentre, thus a
% one-impulse maneuver is carried at the apocentre.
% 5. On final orbit wait to end point.
% Total dv = 2.5498 km/s
% Total time = 43544 s
% Remarks:
% i. The semi-major axis of the first transfer obit has been left as a
% parameter a_t1. The final transfer (after .3) is carried by a
% bitangent approach. The optimal a_t1 has been found in order to minimize
% the total delta-v cost, that corresponds to when the transfer orbit and
% final orbit are tangent at apocentre.
% ii. The suitable point, at which the first tangent burn is carried, is
% chosen such that, after the subsequent plane chane maneuver, the transfer
% orbit and the final orbit are coaxial (same argument of pericentre).
% The calculations of such point are carried numerically in the function
% tangent_transfer_coaxial_orbits.m
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

a_t1 = 14660.6393;   % ~14660.6393 optimal value
[e_t1, th_b1, om_t1, dv_b1] = tangent_transfer_coaxial_orbits([a1 e1 i1 OM1 om1 th1], [a2 e2 i2 OM2 om2 th2], a_t1, mu);
i_t1 = i1;
OM_t1 = OM1;
th_t1 = th_b1 + om1 - om_t1;
dt_1 = timeOfFlight(a1, e1, th1, th_b1, mu);

% Plane change maneuver
[dv_b2, om_t2, th_t2] = changeOrbitalPlane(a_t1, e_t1, i_t1, OM_t1, om_t1, i2, OM2, mu, 1);
a_t2 = a_t1;
e_t2 = e_t1;
i_t2 = i2;
OM_t2 = OM2;
dt_2 = timeOfFlight(a_t1, e_t1, th_t1, th_t2, mu);

% Bitangent transfer between transfer orbit 2 and final orbit
[dv_b3, dv_b4, ~, dt] = bitangent_changeOrbitShape(a_t2, e_t2, om2, a2, e2, om2, mu, 1);
a_t3 = (a2 * ( 1 + e2) + a_t2 * (1 - e_t2)) / 2;
e_t3 = (a2 * ( 1 + e2) - a_t2 * (1 - e_t2)) / 2 / a_t3;
i_t3 = i2;
OM_t3 = OM2;
om_t3 = om2;
dt_3 = timeOfFlight(a_t2, e_t2, th_t2, 0, mu);
dt_4 = pi * sqrt( a_t3 ^ 3 / mu);
dt_5 = timeOfFlight(a2, e2, pi, th2, mu);

dv_tot = abs(dv_b1) + abs(dv_b2) + abs(dv_b3) + abs(dv_b4)
T = dt_1 + dt_2 + dt_3 + dt_4 + dt_5

%% Plot Orbits

% Colors 
c1 = "#48E8C8";
c2 = "#FFBB00";
c3 = "#FF8100";

figure('Name','Alternative Transfer N2', 'NumberTitle', 'Off') 
set(gcf,'color','w');
Terra3d;

[X_t1, Y_t1, Z_t1] = plotOrbit([a_t1, e_t1, i_t1, OM_t1, om_t1, th_t1], mu, 2*pi, pi/1000);
[X_t2, Y_t2, Z_t2] = plotOrbit([a_t2, e_t2, i_t2, OM_t2, om_t2, th_t2], mu, 2 * pi, pi/1000);
plot3(X_t1, Y_t1, Z_t1, '--', 'Color', c1, 'LineWidth', 0.8)
plot3(X_t2, Y_t2, Z_t2,'--', 'Color', c2, 'LineWidth', 0.8)

[X1, Y1, Z1] = plotOrbit([a1, e1, i1, OM1, om1, th1], mu, 2*pi, pi/1000);
[X2, Y2, Z2] = plotOrbit([a2, e2, i2, OM2, om2, th2], mu, 2*pi, pi/1000);
[X_t1, Y_t1, Z_t1] = plotOrbit([a_t1, e_t1, i_t1, OM_t1, om_t1, th_t1], mu, th_t2 - th_t1, pi/1000);
[X_t2, Y_t2, Z_t2] = plotOrbit([a_t2, e_t2, i_t2, OM_t2, om_t2, th_t2], mu, - th_t2, pi/1000);
[X_t3, Y_t3, Z_t3] = plotOrbit([a_t3, e_t3, i_t3, OM_t3, om_t3, 0], mu, pi, pi/1000);

[aX_1, aY_1, aZ_1] = plotOrbit([a1, e1, i1, OM1, om1, 0], mu, pi, pi);
[aX_2, aY_2, aZ_2] = plotOrbit([a2, e2, i2, OM2, om2, 0], mu, pi, pi);
[aX_t1, aY_t1, aZ_t1] = plotOrbit([a_t1, e_t1, i_t1, OM_t1, om_t1, 0], mu, pi, pi);
[aX_t2, aY_t2, aZ_t2] = plotOrbit([a_t2, e_t2, i_t2, OM_t2, om_t2, 0], mu, pi, pi);

plot3(X1, Y1, Z1,'Color','b', 'LineWidth', 1)
plot3(X2, Y2, Z2,'Color','r', 'LineWidth', 1)
plot3(X_t1, Y_t1, Z_t1,'Color', c1, 'LineWidth', 1)
plot3(X_t2, Y_t2, Z_t2,'Color', c2, 'LineWidth', 1)
plot3(X_t3, Y_t3, Z_t3,'Color', c2, 'LineWidth', 1)

plot3(X1(1), Y1(1), Z1(1),'o','Color', 'b','MarkerSize', 6, 'MarkerFaceColor', 'b')
plot3(X2(1), Y2(1), Z2(1),'o','Color', 'r','MarkerSize', 6, 'MarkerFaceColor', 'r')
plot3(X_t1(1), Y_t1(1), Z_t1(1),'x','Color', c1,'MarkerSize', 8, 'LineWidth', 2)
plot3(X_t2(1), Y_t2(1), Z_t2(1),'x','Color', c2,'MarkerSize', 8, 'LineWidth', 2)
plot3(X_t3(end), Y_t3(end), Z_t3(end), 'x','Color', 'r','MarkerSize', 8, 'LineWidth', 2)

plot3(aX_1, aY_1, aZ_1, '-.','Color', 'b','LineWidth', 0.5)
plot3(aX_2, aY_2, aZ_2, '-.','Color', 'r','LineWidth', 0.5)
plot3(aX_t1, aY_t1, aZ_t1, '-.','Color', c1,'LineWidth', 0.5)

grid on
axis equal
box off
title('\fontsize{15}{0}\selectfont Alternative Transfer N2','Interpreter','latex')
subtitle_text = sprintf('dv = %.3f km/s, T = %0.f s', dv_tot, T);
subtitle(subtitle_text, 'Interpreter', 'latex')
xlabel('x [km] - $\gamma$', 'Interpreter', 'latex') 
ylabel('y [km]', 'Interpreter', 'latex') 
zlabel('z [km]', 'Interpreter', 'latex') 
ax = gca;
ax.XRuler.Exponent = 0;
ax.YRuler.Exponent = 0;
ax.ZRuler.Exponent = 0;
