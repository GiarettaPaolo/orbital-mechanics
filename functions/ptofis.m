function [xn, it] = ptofis(x0, phi, nmax, toll)

% ptofis.m - Fixed-point iteration algorithm
%
% PROTOTYPE:
% [xn, it] = ptofis(x0, phi, nmax, toll)
%
% DESCRIPTION:
% Fixed point iterations algorithm to solve phi(x) = x equation
%
% INPUT:
% x0                    [1x1]                   Initial guess                           [-]
% phi                   [1x1->1x1]              Function                                [-]
% nmax                  [1x1]                   Max N. Iterations                       [-]
% toll                  [1x1]                   Tolerance                               [-]
% OUTPUT:
% xn                    [1x1]                   Estimated solution                      [-]
% it                    [1x1]                   Iterations to convergence               [-]
%

err = 1 + toll;
it = 0;
xv = x0;
% Main iterative loop
% New estimate obtained from previous one from x_n+1 = phi(x_n)
% Algorithm halts when maximum number of iterations is reached or the
% absolute difference of consecutive convergents is below the tolerance
while (it < nmax && err > toll)
xn = phi(xv);
err = abs(xn - xv);
it = it + 1;
xv = xn;
end