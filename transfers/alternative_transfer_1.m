% Alternative transfer N°1
% The transfer consists of the following maneuvers:
% 1. On initial orbit wait to maneuver point.
% 2. Tangent burn at pericentre to get to desired transfer orbit.
% 3. Plane change maneuver at minimum delta-v cost point.
% 4. Secant transfer to final orbit at minimum delta_v cost point.
% 5. On final orbit wait to end point.
% Total dv = 2.4763 km/s
% Total time = 23756 s
% Remarks:
% i. The semi major-axis of the first transfer orbit has been left as a
% parameter (a_t1). For the proposed solution the semi-major axis has been
% optimized for minimum total delta-v.
% ii. The semi-major axis of the first transfer orbit must be large enough
% to guarantee an intersection with the final orbit for the secant
% transfer.
% iii. The optimal transfer procedure is obtained when the two secant
% orbits are close in being tangent, but not perfectly tangent.
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

% Pericentre tangent burn
a_t1 = 14019.27325714;  % Optimal value 14019.27325714
e_t1 = 1 - a1 * (1 - e1) / a_t1;
i_t1 = i1;
OM_t1 = OM1;
om_t1 = om1;
% delta-v calculated from velocity difference as burn is tangent
dv_b1 = sqrt(mu / a_t1 / (1 - e_t1^2)) * (1 + e_t1) - sqrt(mu / a1 / (1 - e1^2)) * (1 + e1); 
th_t1 = 0;
dt_1 = timeOfFlight(a1, e1, th1, th_t1, mu);

% Plane Change
[dv_b2, om_t2, th_t2] = changeOrbitalPlane(a_t1, e_t1, i_t1, OM_t1, om_t1, i2, OM2, mu, 1);
a_t2 = a_t1;
e_t2 = e_t1;
i_t2 = i2;
OM_t2 = OM2;
dt_2 = timeOfFlight(a_t1, e_t1, th_t1, th_t2, mu);

% Secant solution between transfer orbit 2 and final orbit
% The equation can be written as a linear trigonometric equation that is
% solved below by passing in polar form.
p2 = a2 * (1 - e2^2);
p_t2 = a_t2 * (1 - e_t2^2);

%C = A * cos(th) + B * sin(th) = G cos(th + phi)
C = p2 - p_t2;
A = (e2 * p_t2 * cos(om2 - om_t2) - e_t2 * p2);
B = e2 * p_t2 * sin(om2 - om_t2);

G = sqrt(A^2 + B^2);
phi = atan2(-B/G, A/G);

% Check that an intersection exists
if abs(C) > abs(G)
    error('No intersection between orbits')
end

th_b3 = mod(acos(C/G) - phi, 2 * pi); % <- put a - sign infront of the acos for the other (sub-optimal) secant solution
th_t3 = mod(th_b3 + om_t2 - om2, 2 * pi);

% Time of flight to get to intersection point.
dt_3 = timeOfFlight(a_t1, e_t1, th_t2, th_b3, mu);

% Delta-v calculation from radial and tangential velocities before and
% after the impulsive maneuver.
vr1 = sqrt(mu / p_t2) * e_t2 * sin(th_b3);
vt1 = sqrt(mu / p_t2) * (1 + e_t2 * cos(th_b3));
vr2 = sqrt(mu / p2) * e2 * sin(th_t3);
vt2 = sqrt(mu / p2) * (1 + e2 * cos(th_t3));
dv_b3 = sqrt((vr2 - vr1) ^ 2 + (vt2 - vt1) ^ 2);

% Time of flight to get to final point
dt_4 = timeOfFlight(a2, e2, th_t3, th2, mu);

T = dt_1 + dt_2 + dt_3 + dt_4
dv_tot = abs(dv_b1) + abs(dv_b2) + abs(dv_b3)

%% Plot Orbits

% Colors 
c1 = "#48E8C8";
c2 = "#FFBB00";
c3 = "#FF8100";

figure('Name','Alternative Transfer N1', 'NumberTitle', 'Off') 
set(gcf,'color','w');
Terra3d;

[X_t1, Y_t1, Z_t1] = plotOrbit([a_t1, e_t1, i_t1, OM_t1, om_t1, th_t1], mu, 2*pi, pi/1000);
[X_t2, Y_t2, Z_t2] = plotOrbit([a_t2, e_t2, i_t2, OM_t2, om_t2, th_t2], mu, 2 * pi, pi/1000);
plot3(X_t1, Y_t1, Z_t1, '--', 'Color', c1, 'LineWidth', 0.8)
plot3(X_t2, Y_t2, Z_t2,'--', 'Color', c2, 'LineWidth', 0.8)

[X1, Y1, Z1] = plotOrbit([a1, e1, i1, OM1, om1, th1], mu, 2*pi, pi/1000);
[X2, Y2, Z2] = plotOrbit([a2, e2, i2, OM2, om2, th2], mu, 2*pi, pi/1000);
[X_t1, Y_t1, Z_t1] = plotOrbit([a_t1, e_t1, i_t1, OM_t1, om_t1, th_t1], mu, th_t2 - th_t1, pi/1000);
[X_t2, Y_t2, Z_t2] = plotOrbit([a_t2, e_t2, i_t2, OM_t2, om_t2, th_t2], mu, th_b3 - th_t2, pi/1000);

[aX_1, aY_1, aZ_1] = plotOrbit([a1, e1, i1, OM1, om1, 0], mu, pi, pi);
[aX_2, aY_2, aZ_2] = plotOrbit([a2, e2, i2, OM2, om2, 0], mu, pi, pi);
[aX_t1, aY_t1, aZ_t1] = plotOrbit([a_t1, e_t1, i_t1, OM_t1, om_t1, 0], mu, pi, pi);
[aX_t2, aY_t2, aZ_t2] = plotOrbit([a_t2, e_t2, i_t2, OM_t2, om_t2, 0], mu, pi, pi);

plot3(X1, Y1, Z1,'Color','b', 'LineWidth', 1)
plot3(X2, Y2, Z2,'Color','r', 'LineWidth', 1)
plot3(X_t1, Y_t1, Z_t1,'Color', c1, 'LineWidth', 1)
plot3(X_t2, Y_t2, Z_t2,'Color', c2, 'LineWidth', 1)

plot3(X1(1), Y1(1), Z1(1),'o','Color', 'b','MarkerSize', 6, 'MarkerFaceColor', 'b')
plot3(X2(1), Y2(1), Z2(1),'o','Color', 'r','MarkerSize', 6, 'MarkerFaceColor', 'r')
plot3(X_t1(1), Y_t1(1), Z_t1(1),'x','Color', c1,'MarkerSize', 8, 'LineWidth', 2)
plot3(X_t2(1), Y_t2(1), Z_t2(1),'x','Color', c2,'MarkerSize', 8, 'LineWidth', 2)
plot3(X_t2(end), Y_t2(end), Z_t2(end), 'x','Color', 'r','MarkerSize', 8, 'LineWidth', 2)

plot3(aX_1, aY_1, aZ_1, '-.','Color', 'b','LineWidth', 0.5)
plot3(aX_2, aY_2, aZ_2, '-.','Color', 'r','LineWidth', 0.5)
plot3(aX_t1, aY_t1, aZ_t1, '-.','Color', c1,'LineWidth', 0.5)

grid on
axis equal
box off
title('\fontsize{15}{0}\selectfont Alternative Transfer N1','Interpreter','latex')
subtitle_text = sprintf('dv = %.3f km/s, T = %0.f s', dv_tot, T);
subtitle(subtitle_text, 'Interpreter', 'latex')
xlabel('x [km] - $\gamma$', 'Interpreter', 'latex') 
ylabel('y [km]', 'Interpreter', 'latex') 
zlabel('z [km]', 'Interpreter', 'latex') 
ax = gca;
ax.XRuler.Exponent = 0;
ax.YRuler.Exponent = 0;
ax.ZRuler.Exponent = 0;
