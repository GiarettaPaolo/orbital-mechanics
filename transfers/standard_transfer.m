% Standard orbital transfer
% The standard orbital transfer consists of the following maneuvers:
% 1. On initial orbit wait to maneuver point.
% 2. Change plane maneuver to get coplanar with final orbit.
% 3. Change argument of pericentre to get coaxial with the final orbit.
% 4. Pericentre - Apocentre bitanangent transfer to final orbit.
% 5. On final orbit wait to end point.
% Total dv = 3.2815 km/s (4 impulsive burns)
% Total time = 34274 s
% Remarks:
% i. The change of argument of pericentre maneuver has the same delta-v cost
% independently of the choice of maneuver point, thus the first point
% encountered by the satellite in in its trajectory has been used.
% ii. The change of plane maneuver might have different delta-v costs
% depending on the maneuver point, however, in this case the first
% intersection point is also the one with the minimum cost. 
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

% Plane Change maneuver 
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

% Time transfer to get to maneuver point 3 (pericente)
th_b3 = 0;
dt_3 = timeOfFlight(a_t2, e_t2, th_t2, th_b3);

% Bitangent transfer pericenter - apocenter
[dv_b3, dv_b4, th_b4, dt_4] = bitangent_changeOrbitShape(a_t2, e_t2, om_t2, a2, e2, om2, mu, 1);
a_t3 = (a_t2 * (1 - e_t2) + a2 * (1 + e2)) / 2;
e_t3 = (a2 * (1 + e2) - a_t2 * (1 - e_t2)) / 2 / a_t3;
OM_t3 = OM2;
om_t3 = om2;
i_t3 = i2;


% Time transfer to get to final point
dt_5 = timeOfFlight(a2, e2, th_b4, th2, mu);

% Total time
T = dt_1 + dt_2 + dt_3 + dt_4 + dt_5
% Total delta V
dv_tot = abs(dv_b1) + abs(dv_b2) + abs(dv_b3) + abs(dv_b4)


%% Plot Orbits

% Colors 
c1 = "#48E8C8";
c2 = "#FFBB00";
c3 = "#FF8100";

% First figure
figure('Name','Standard Transfer', 'NumberTitle', 'Off') 
set(gcf,'color','w');
Terra3d;

[X_t1, Y_t1, Z_t1] = plotOrbit([a_t1, e_t1, i_t1, OM_t1, om_t1, th_t1], mu, 2*pi, pi/1000);
[X_t2, Y_t2, Z_t2] = plotOrbit([a_t2, e_t2, i_t2, OM_t2, om_t2, th_t2], mu, 2 * pi, pi/1000);
[X_t3, Y_t3, Z_t3] = plotOrbit([a_t3, e_t3, i_t3, OM_t3, om_t3, 0], mu, 2*pi , pi/1000);
plot3(X_t1, Y_t1, Z_t1, '--', 'Color', c1, 'LineWidth', 0.8)
plot3(X_t2, Y_t2, Z_t2,'--', 'Color', c2, 'LineWidth', 0.8)
plot3(X_t3, Y_t3, Z_t3,'--', 'Color', c3, 'LineWidth', 0.8)


[X1, Y1, Z1] = plotOrbit([a1, e1, i1, OM1, om1, th1], mu, 2*pi, pi/1000);
[X2, Y2, Z2] = plotOrbit([a2, e2, i2, OM2, om2, th2], mu, 2*pi, pi/1000);
[X_t1, Y_t1, Z_t1] = plotOrbit([a_t1, e_t1, i_t1, OM_t1, om_t1, th_t1], mu, th_b2 - th_t1, pi/1000);
[X_t2, Y_t2, Z_t2] = plotOrbit([a_t2, e_t2, i_t2, OM_t2, om_t2, th_t2], mu, th_b3 - th_t2, pi/1000);
[X_t3, Y_t3, Z_t3] = plotOrbit([a_t3, e_t3, i_t3, OM_t3, om_t3, 0], mu, th_b4, pi/1000);

[aX_1, aY_1, aZ_1] = plotOrbit([a1, e1, i1, OM1, om1, 0], mu, pi, pi);
[aX_2, aY_2, aZ_2] = plotOrbit([a2, e2, i2, OM2, om2, 0], mu, pi, pi);
[aX_t1, aY_t1, aZ_t1] = plotOrbit([a_t1, e_t1, i_t1, OM_t1, om_t1, 0], mu, pi, pi);
[aX_t2, aY_t2, aZ_t2] = plotOrbit([a_t2, e_t2, i_t2, OM_t2, om_t2, 0], mu, pi, pi);

plot3(X1, Y1, Z1,'Color','b', 'LineWidth', 1)
plot3(X2, Y2, Z2,'Color','r', 'LineWidth', 1)
plot3(X_t1, Y_t1, Z_t1,'Color', c1, 'LineWidth', 1)
plot3(X_t2, Y_t2, Z_t2,'Color', c2, 'LineWidth', 1)
plot3(X_t3, Y_t3, Z_t3,'Color', c3, 'LineWidth', 1)

plot3(X1(1), Y1(1), Z1(1),'o','Color', 'b','MarkerSize', 6, 'MarkerFaceColor', 'b')
plot3(X2(1), Y2(1), Z2(1),'o','Color', 'r','MarkerSize', 6, 'MarkerFaceColor', 'r')
plot3(X_t1(1), Y_t1(1), Z_t1(1),'x','Color', c1,'MarkerSize', 8, 'LineWidth', 2)
plot3(X_t2(1), Y_t2(1), Z_t2(1),'x','Color', c2,'MarkerSize', 8, 'LineWidth', 2)
plot3(X_t3(1), Y_t3(1), Z_t3(1),'x','Color', c3,'MarkerSize', 8, 'LineWidth', 2)
plot3(X_t3(end), Y_t3(end), Z_t3(end), 'x','Color', 'r','MarkerSize', 8, 'LineWidth', 2)

plot3(aX_1, aY_1, aZ_1, '-.','Color', 'b','LineWidth', 0.5)
plot3(aX_2, aY_2, aZ_2, '-.','Color', 'r','LineWidth', 0.5)
plot3(aX_t1, aY_t1, aZ_t1, '-.','Color', c1,'LineWidth', 0.5)

grid on
axis equal
box off
title('\fontsize{15}{0}\selectfont Standard Transfer','Interpreter','latex')
subtitle_text = sprintf('dv = %.3f km/s, T = %0.f s', dv_tot, T);
subtitle(subtitle_text, 'Interpreter', 'latex')
xlabel('x [km] - $\gamma$', 'Interpreter', 'latex') 
ylabel('y [km]', 'Interpreter', 'latex') 
zlabel('z [km]', 'Interpreter', 'latex') 
ax = gca;
ax.XRuler.Exponent = 0;
ax.YRuler.Exponent = 0;
ax.ZRuler.Exponent = 0;

% Second figure
figure('Name','Standard Transfer - Coplanar Transfers', 'NumberTitle', 'Off') 
set(gcf,'color','w');
hold on

[X_t1, Y_t1, ~] = plotOrbit([a_t1, e_t1, 0, 0, om_t1 - om2, th_t1], mu, 2*pi, pi/1000);
[X_t2, Y_t2, ~] = plotOrbit([a_t2, e_t2, 0, 0, om_t2 - om2, th_t2], mu, 2 * pi, pi/1000);
[X_t3, Y_t3, ~] = plotOrbit([a_t3, e_t3, 0, 0, om_t3 - om2, 0], mu, 2*pi , pi/1000);
plot(X_t1, Y_t1, '--', 'Color', c1, 'LineWidth', 0.8)
plot(X_t2, Y_t2, '--', 'Color', c2, 'LineWidth', 0.8)
plot(X_t3, Y_t3, '--', 'Color', c3, 'LineWidth', 0.8)

[X2, Y2, ~] = plotOrbit([a2, e2, 0, 0, 0, th2], mu, 2*pi, pi/1000);
[X_t1, Y_t1, ~] = plotOrbit([a_t1, e_t1, 0, 0, om_t1-om2, th_t1], mu, th_b2 - th_t1, pi/1000);
[X_t2, Y_t2, ~] = plotOrbit([a_t2, e_t2, 0, 0, om_t2-om2, th_t2], mu, th_b3 - th_t2, pi/1000);
[X_t3, Y_t3, ~] = plotOrbit([a_t3, e_t3, 0, 0, om_t3-om2, 0], mu, th_b4, pi/1000);
plot(X2, Y2, '-', 'Color', 'r', 'LineWidth', 1.2)
plot(X_t1, Y_t1, '-', 'Color', c1, 'LineWidth', 1.2)
plot(X_t2, Y_t2, '-', 'Color', c2, 'LineWidth', 1.2)
plot(X_t3, Y_t3, '-', 'Color', c3, 'LineWidth', 1.2)

[aX_2, aY_2, ~] = plotOrbit([a2, e2, 0, 0, 0, 0], mu, pi, pi);
[aX_t1, aY_t1, ~] = plotOrbit([a_t1, e_t1, 0, 0, om_t1-om2, 0], mu, pi, pi);
plot(aX_2, aY_2, '-.','Color', 'r','LineWidth', 0.5)
plot(aX_t1, aY_t1, '-.','Color', c1,'LineWidth', 0.5)

plot(X2(1), Y2(1),'o','Color', 'r','MarkerSize', 6, 'MarkerFaceColor', 'r')
plot(X_t1(1), Y_t1(1),'x','Color', c1,'MarkerSize', 8, 'LineWidth', 2)
plot(X_t2(1), Y_t2(1),'x','Color', c2,'MarkerSize', 8, 'LineWidth', 2)
plot(X_t3(1), Y_t3(1),'x','Color', c3,'MarkerSize', 8, 'LineWidth', 2)
plot(X_t3(end), Y_t3(end), 'x','Color', 'r','MarkerSize', 8, 'LineWidth', 2)

grid on 
axis equal
box off
title('\fontsize{15}{0}\selectfont Standard Transfer - Coplanar Transfers','Interpreter','latex')
subtitle_text = sprintf('Perifocal reference frame of final orbit');
subtitle(subtitle_text, 'Interpreter', 'latex')
xlabel('$r_{\hat{e}}$ [km]' , 'Interpreter', 'latex') 
ylabel('$r_{\hat{p}}$ [km]', 'Interpreter', 'latex')  
ax = gca;
ax.XRuler.Exponent = 0;
ax.YRuler.Exponent = 0;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';

