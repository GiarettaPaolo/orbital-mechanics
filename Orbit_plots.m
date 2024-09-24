clear all
close all
clc

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

% Plot initial / final orbit
figure('Name','Assigned Orbit', 'NumberTitle', 'Off') 
set(gcf,'color','w');
hold on

[X1, Y1, Z1] = plotOrbit([a1, e1, i1, OM1, om1, th1], mu, 2*pi, pi/1000);
[X2, Y2, Z2] = plotOrbit([a2, e2, i2, OM2, om2, th2], mu, 2*pi, pi/1000);

[aX_1, aY_1, aZ_1] = plotOrbit([a1, e1, i1, OM1, om1, 0], mu, pi, pi);
[aX_2, aY_2, aZ_2] = plotOrbit([a2, e2, i2, OM2, om2, 0], mu, pi, pi);

Terra3d;
h(1) = plot3(X1, Y1, Z1,'Color','b', 'LineWidth', 1.5, 'DisplayName', 'Initial Orbit');
h(2) = plot3(X2, Y2, Z2,'Color','r', 'LineWidth', 1.5, 'DisplayName', 'Final Orbit');
plot3(X1(1), Y1(1), Z1(1),'o','Color', 'b','MarkerSize', 10, 'MarkerFaceColor', 'b')
plot3(X2(1), Y2(1), Z2(1),'o','Color', 'r','MarkerSize', 10, 'MarkerFaceColor', 'r')
plot3(aX_1, aY_1, aZ_1, '-.','Color', 'b','LineWidth', 0.5)
plot3(aX_2, aY_2, aZ_2, '-.','Color', 'r','LineWidth', 0.5)

grid on 
axis equal
box off
title('\fontsize{15}{0}\selectfont Assigned orbits','Interpreter','latex')
xlabel('x [km] - $\gamma$', 'Interpreter', 'latex') 
ylabel('y [km]', 'Interpreter', 'latex') 
zlabel('z [km]', 'Interpreter', 'latex') 
legend(h(1:2))
ax = gca;
ax = gca;
ax.XRuler.Exponent = 0;
ax.YRuler.Exponent = 0;
ax.ZRuler.Exponent = 0;
view(0, 0)

%%

figure('Name','First Orbit - Perifocal View', 'NumberTitle', 'Off') 
set(gcf,'color','w');
hold on

[X1, Y1, ~] = plotOrbit([a1, e1, 0, 0, 0, th1], mu, 2*pi, pi/1000);
plot(X1, Y1, '-', 'Color', 'b', 'LineWidth', 1.2)
plot(X1(1), Y1(1),'o','Color', 'b','MarkerSize', 6, 'MarkerFaceColor', 'b')


grid on 
axis equal
box off
title('\fontsize{15}{0}\selectfont Initial Orbit - Perifocal View','Interpreter','latex')
xlabel('$r_{\hat{e}}$ [km]' , 'Interpreter', 'latex') 
ylabel('$r_{\hat{p}}$ [km]', 'Interpreter', 'latex')  
ax = gca;
ax.XRuler.Exponent = 0;
ax.YRuler.Exponent = 0;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';

%%

figure('Name','Final Orbit - Perifocal View', 'NumberTitle', 'Off') 
set(gcf,'color','w');
hold on

[X2, Y2, ~] = plotOrbit([a2, e2, 0, 0, 0, th2], mu, 2*pi, pi/1000);
plot(X2, Y2, '-', 'Color', 'r', 'LineWidth', 1.2)
plot(X2(1), Y2(1),'o','Color', 'r','MarkerSize', 6, 'MarkerFaceColor', 'r')


grid on 
axis equal
box off
title('\fontsize{15}{0}\selectfont Final Orbit - Perifocal View','Interpreter','latex')
xlabel('$r_{\hat{e}}$ [km]' , 'Interpreter', 'latex') 
ylabel('$r_{\hat{p}}$ [km]', 'Interpreter', 'latex')  
ax = gca;
ax.XRuler.Exponent = 0;
ax.YRuler.Exponent = 0;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';

