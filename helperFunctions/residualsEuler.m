function [continuityResidual, xMomentumResidual, yMomentumResidual] = residualsEuler(U1, U2, dx, dy)
% U1 - discretised flow field variables (rho, u, v, P) elements 1-4 of 3rd dimension of iteration i
% U2 - discretised flow field variables (rho, u, v, P) elements 1-4 of 3rd dimension of iteration i
% dx, dy - grid spacing in x, y ----- THESE MUST BE CONSTANTS -----

[rho1, u1, v1, p1] = deal(U1(:,:,1), U1(:,:,2), U1(:,:,3), U1(:,:,4));
[rho2, u2, v2, p2] = deal(U2(:,:,1), U1(:,:,2), U1(:,:,3), U1(:,:,4));

% continuity eqn.
% du/dx + dv/dy = 0 
continuityResidual = max((rho2 .* u2 - rho1 .* u1)/dx + (rho2 .* v2 - rho1 .* v1)/dy,[], "all"); 

% x-momentum eqn.
% (d rho u^2)/dx + (d rho uv)/dy + dp/dx = 0
d_rho_u2 = ((rho2 .* u2 .^ 2) - (rho1 .* u1 .^ 2));
d_rho_uv = ((rho2 .* u2 .* v2) - (rho1 .* u1 .* v1));
dp =  p2 - p1;
xMomentumResidual = max(d_rho_u2/dx + d_rho_uv/dy + dp/dx,[], "all");

% y-momentum eqn.
% (d rho uv)/dx + (d rho v^2)/dy + dp/dy = 0
d_rho_v2 = ((rho2 .* v2 .^ 2) - (rho1 .* v1 .^ 2));
yMomentumResidual =  max(d_rho_uv/dx + d_rho_v2/dy + dp/dy,[], "all");

% EoS (Ideal Gas)
% P - rho * R * T = 0

idealGasResidual = (p2 - p1);

end