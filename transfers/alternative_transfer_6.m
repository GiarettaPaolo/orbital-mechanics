% Alternative transfer N°6 - Optimal two-impulse maneuver
% The transfer consists of the following maneuvers:
% 1. On initial orbit wait to maneuver point.
% 2. First impulsive maneuver.
% 3. Second impulsive maneuver.
% 4. On final orbit wait to end point.
% Total dv = 2.0735 km/S
% Total time = 23631 s
% Remarks:
% i. Between the infinitude of transfer orbits, those that would result
% in the satellite intersecting the Earth's surface were excluded.
% ii. The proposed solution is the orbit that minimizes the total delta-v
% cost of the transfer by using the function direct_optimal_transfer.m and
% finding the first and second burn points that minimize the total cost.
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

% Main optimization loop (find orbit with minimum delta_v)
tol = 1e-9;            % tolerance
n = 20;                % number of sample points for iteration
step = 1 + tol;     
c1 = 0;                % initial interval spanned for th_b1
b1 = 2*pi;             % is [c1, b1]
c2 = 0;                % initial interval spanned for th_b2
b2 = 2*pi;             % is [c2, b2]
dv_min = inf;
th1_min = 0;          
th2_min = 0;

while step > tol
    
    % Step size for every iteration
    step = (b1 - c1) / n;          
    
    for th_b1 = c1:step:b1
        for th_b2 = c2:step:b2
            
            % Delta-v cost of the direct transfer between points
            [a_t, e_t, i_t, OM_t, om_t, th_t1, th_t2, dv_tot] = direct_optimal_transfer([a1, e1, i1, OM1, om1, th_b1], ...
                [a2, e2, i2, OM2, om2, th_b2]);
            
            % Local minimum of values evalued at sample points
            if dv_tot < dv_min
                dv_min = dv_tot;
                th1_min = th_b1;
                th2_min = th_b2;
            end
            
        end
    end
    
    % New interval to be sampled for the next iteration
    c1 = th1_min - step;
    b1 = th1_min + step;
    c2 = th2_min - step;
    b2 = th2_min + step;
    
end

% Time of flight
dt_1 = timeOfFlight(a1, e1, th1, th_b1, mu);
dt_2 = timeOfFlight(a_t, e_t, th_t1, th_t2, mu);
dt_3 = timeOfFlight(a2, e2, th_b2, th2, mu);

% Delta-v calculated as norm of velocity vectors difference
[r_t1, v_t1] = kep2car(a_t, e_t, i_t, OM_t, om_t, th_t1, mu);
[r_t2, v_t2] = kep2car(a_t, e_t, i_t, OM_t, om_t, th_t2, mu);
[r_b1, v_b1] = kep2car(a1, e1, i1, OM1, om1, th_b1, mu);
[r_b2, v_b2] = kep2car(a2, e2, i2, OM2, om2, th_b2, mu);
dv_b1 = norm(v_t1 - v_b1);
dv_b2 = norm(v_t2 - v_b2);

T = dt_1 + dt_2 + dt_3
dv_tot

%% Plot Orbits

% Colors 
c1 = "#48E8C8";
c2 = "#FFBB00";
c3 = "#FF8100";

figure('Name','Alternative Transfer N6', 'NumberTitle', 'Off')
set(gcf,'color','w');
Terra3d;

[X_t1, Y_t1, Z_t1] = plotOrbit([a_t, e_t, i_t, OM_t, om_t, th_t1], mu, 2*pi, pi/1000);
plot3(X_t1, Y_t1, Z_t1, '--', 'Color', c3, 'LineWidth', 0.8)

[X1, Y1, Z1] = plotOrbit([a1, e1, i1, OM1, om1, th1], mu, 2*pi, pi/1000);
[X2, Y2, Z2] = plotOrbit([a2, e2, i2, OM2, om2, th2], mu, 2*pi, pi/1000);
[X_t, Y_t, Z_t] = plotOrbit([a_t, e_t, i_t, OM_t, om_t, th_t1], mu, th_t2 - th_t1, pi/1000);

[aX_1, aY_1, aZ_1] = plotOrbit([a1, e1, i1, OM1, om1, 0], mu, pi, pi);
[aX_2, aY_2, aZ_2] = plotOrbit([a2, e2, i2, OM2, om2, 0], mu, pi, pi);
[aX_t, aY_t, aZ_t] = plotOrbit([a_t, e_t, i_t, OM_t, om_t 0], mu, pi, pi);


plot3(X1, Y1, Z1,'Color','b', 'LineWidth', 1)
plot3(X2, Y2, Z2,'Color','r', 'LineWidth', 1)
plot3(X_t, Y_t, Z_t,'Color', c2, 'LineWidth', 1)

plot3(X1(1), Y1(1), Z1(1),'o','Color', 'b','MarkerSize', 6, 'MarkerFaceColor', 'b')
plot3(X2(1), Y2(1), Z2(1),'o','Color', 'r','MarkerSize', 6, 'MarkerFaceColor', 'r')
plot3(X_t(1), Y_t(1), Z_t(1),'x','Color', c2,'MarkerSize', 8, 'LineWidth', 2)
plot3(X_t(end), Y_t(end), Z_t(end), 'x','Color', 'r','MarkerSize', 8, 'LineWidth', 2)

plot3(aX_1, aY_1, aZ_1, '-.','Color', 'b','LineWidth', 0.5)
plot3(aX_2, aY_2, aZ_2, '-.','Color', 'r','LineWidth', 0.5)
plot3(aX_t, aY_t, aZ_t, '-.','Color', c2,'LineWidth', 0.5)

grid on
axis equal
box off
title('\fontsize{15}{0}\selectfont Alternative Transfer N6','Interpreter','latex')
subtitle_text = sprintf('dv = %.3f km/s, T = %0.f s', dv_tot, T);
subtitle(subtitle_text, 'Interpreter', 'latex')
xlabel('x [km] - $\gamma$', 'Interpreter', 'latex') 
ylabel('y [km]', 'Interpreter', 'latex') 
zlabel('z [km]', 'Interpreter', 'latex') 
ax = gca;
ax.XRuler.Exponent = 0;
ax.YRuler.Exponent = 0;
ax.ZRuler.Exponent = 0;