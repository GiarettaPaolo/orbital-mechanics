% extra transfer 5
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

% Main optimization loop (find orbit with minimum delta_v)
tol = 1e-12;
n = 40;
step = 1 + tol;
c1 = 0;
b1 = 2*pi;
c2 = 10000;
b2 = 100000;
dv_min = inf;
th_min = 0;
a_t1= 0;

while step > tol
    
    step = (b1 - c1) / n;
    
    for th_t1 = c1:step:b1
        for a_t1 = c2:step:b2


        % Tangent burn
        % Optimal value 14019.27325714 -> dv=2.47632716174967
        
        [e_t1, om_t1, dv_b1] = tangent_burn(a1, e1, om1, th_t1, a_t1, mu);
        i_t1 = i1;
        OM_t1 = OM1;

        % Plane Change
        [dv_b2, om_t2, th_t2] = changeOrbitalPlane(a_t1, e_t1, i_t1, OM_t1, om_t1, i2, OM2, mu, 1);
        a_t2 = a_t1;
        e_t2 = e_t1;
        i_t2 = i2;
        OM_t2 = OM2;

        % Secant solution
        p2 = a2 * (1 - e2^2);
        p_t2 = a_t2 * (1 - e_t2^2);

        %C = A * cos(th) + B * sin(th) = G cos(th + phi)
        C = p2 - p_t2;
        A = (e2 * p_t2 * cos(om2 - om_t2) - e_t2 * p2);
        B = e2 * p_t2 * sin(om2 - om_t2);

        G = sqrt(A^2 + B^2);
        phi = atan2(-B/G, A/G);
        
        
        if abs(C) > abs(G)
            continue
        end
        
        th_b3 = mod(acos(C/G) - phi, 2 * pi); % <- put a - sign infront of the acos for the other (sub_optimal) secant solution
        th_t3 = mod(th_b3 + om_t2 - om2, 2 * pi);

        % Delta-v calculation
        vr1 = sqrt(mu / p_t2) * e_t2 * sin(th_b3);
        vt1 = sqrt(mu / p_t2) * (1 + e_t2 * cos(th_b3));
        vr2 = sqrt(mu / p2) * e2 * sin(th_t3);
        vt2 = sqrt(mu / p2) * (1 + e2 * cos(th_t3));
        dv_b3 = sqrt((vr2 - vr1) ^ 2 + (vt2 - vt1) ^ 2);

        dv_tot = abs(dv_b1) + abs(dv_b2) + abs(dv_b3);

        if dv_tot < dv_min
                    dv_min = dv_tot;
                    th_min = th_t1;
                    a_min = a_t1;
         end

        end
     end
    
    c1 = th_min - step;
    b1 = th_min + step;
    c2 = a_min - step;
    b2 = a_min + step;
    
end

%% Plot Orbits

[X1, Y1, Z1] = plotOrbit([a1, e1, i1, OM1, om1, th1], mu, 2*pi, pi/1000);
[X2, Y2, Z2] = plotOrbit([a2, e2, i2, OM2, om2, th2], mu, 2*pi, pi/1000);
[X_t1, Y_t1, Z_t1] = plotOrbit([a_t1, e_t1, i_t1, OM_t1, om_t1, (th_t1 + om1 - om_t1)], mu, th_t2 - th_t1 - om1 + om_t1, pi/1000);
[X_t2, Y_t2, Z_t2] = plotOrbit([a_t2, e_t2, i_t2, OM_t2, om_t2, th_t2], mu, th_b3 - th_t2, pi/1000);

[aX_1, aY_1, aZ_1] = plotOrbit([a1, e1, i1, OM1, om1, 0], mu, pi, pi);
[aX_2, aY_2, aZ_2] = plotOrbit([a2, e2, i2, OM2, om2, 0], mu, pi, pi);
[aX_t1, aY_t1, aZ_t1] = plotOrbit([a_t1, e_t1, i_t1, OM_t1, om_t1, 0], mu, pi, pi);
[aX_t2, aY_t2, aZ_t2] = plotOrbit([a_t2, e_t2, i_t2, OM_t2, om_t2, 0], mu, pi, pi);

Terra3d;
plot3(X1, Y1, Z1, 'b', 'LineWidth', 1)
plot3(X2, Y2, Z2, 'r', 'LineWidth', 1)
plot3(X_t1, Y_t1, Z_t1, 'c', 'LineWidth', 1)
plot3(X_t2, Y_t2, Z_t2, 'k', 'LineWidth', 1)

plot3(X1(1), Y1(1), Z1(1), 'ob', 'LineWidth', 2)
plot3(X2(1), Y2(1), Z2(1), 'or', 'LineWidth', 2)
plot3(X_t1(1), Y_t1(1), Z_t1(1), 'xc', 'LineWidth', 2)
plot3(X_t2(1), Y_t2(1), Z_t2(1), 'xk', 'LineWidth', 2)
plot3(X_t2(end), Y_t2(end), Z_t2(end), 'xr', 'LineWidth', 2)

plot3(aX_1, aY_1, aZ_1, '-.b', 'LineWidth', 0.5)
plot3(aX_2, aY_2, aZ_2, '-.r', 'LineWidth', 0.5)
plot3(aX_t1, aY_t1, aZ_t1, '-.c', 'LineWidth', 0.5)

grid on